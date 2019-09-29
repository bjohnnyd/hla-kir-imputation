from os import path
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
from snakemake.utils import format as snakeformat

HTTP = HTTPRemoteProvider()

URL_DICT = {
    "HG18" : "http://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz", "HG17" : "https://hgdownload.soe.ucsc.edu/goldenPath/hg17/liftOver/hg17ToHg19.over.chain.gz",
    "HG16" : "https://hgdownload.soe.ucsc.edu/goldenPath/hg16/liftOver/hg16ToHg19.over.chain.gz",
}

rule get_overchain: 
    input: lambda wc: HTTP.remote(URL_DICT[wc.reference.upper()], keep_local = True)
    output: "input/meta/liftover/{reference}ToHg19.over.chain.gz"
    shell: "mv {input} {output}"

rule generate_ucsc_bed:
    input:
        plink = lambda wc: config['project'][wc.project]['liftover']['plink'],
    output:
        plink = "output/{project}/liftover/01_ucsc_bed/{project}.{reference}.map",
        bed = "output/{project}/liftover/01_ucsc_bed/{project}.{reference}.ucsc.bed",
    params:
        basein = lambda wildcards,input: path.splitext(input.plink)[0], baseout =  "output/{project}/liftover/01_ucsc_bed/{project}.{reference}"
    conda: "../envs/liftover.yaml"
    log: "logs/liftover/{project}/{project}.{reference}.01_generate_ucsc_bed.log",
    shell: 
        r"""
        plink --bfile {params.basein} --allow-no-sex --real-ref-alleles --recode --out {params.baseout} &> {log}
        awk 'BEGIN{{OFS="\t"}} {{print "chr"$1,$4-1,$4,$2,$3}}' {output.plink} > {output.bed}
        """
rule liftOver:
    input: 
    	bed = rules.generate_ucsc_bed.output.bed,
	chain = "input/meta/liftover/{reference}ToHg19.over.chain.gz",
    output: 
        lifted = "output/{project}/liftover/02_liftover/{project}.{reference}ToHg19.ucsc.bed",
        unlifted = "output/{project}/liftover/02_liftover/{project}.{reference}ToHg19.unlifted.bed",
    conda: "../envs/liftover.yaml"
    log: "logs/liftover/{project}/{project}.{reference}.02_liftOver.log"
    shell: "liftOver {input.bed} {input.chain} {output.lifted} {output.unlifted} &> {log}"

rule create_plink_lifted:
    input: 
        lifted = rules.liftOver.output.lifted,
        unlifted = rules.liftOver.output.unlifted
    output: 
        lifted_changefile = temp("output/{project}/liftover/03_hg19_plink/{project}.{reference}ToHg19.lifted.list"),
	skip_snps = "output/{project}/liftover/03_hg19_plink/{project}.{reference}.skip_snps.list",
        plink = "output/{project}/liftover/03_hg19_plink/{project}.{reference}ToHg19.plink.bed"
    params:
        basein = "output/{project}/liftover/01_ucsc_bed/{project}.{reference}",
        baseout =  "output/{project}/liftover/03_hg19_plink/{project}.{reference}ToHg19.plink"
    conda: "../envs/liftover.yaml"
    log: "logs/liftover/{project}/{project}.{reference}.03_create_plink_lifted.log"
    shell: 
        r"""
        awk '{{print $4,$3}}' {input.lifted} > {output.lifted_changefile}
        awk '{{print $4}}' {input.unlifted} > {output.skip_snps}
        plink --file {params.basein} --allow-no-sex --real-ref-alleles --make-bed --update-map {output.lifted_changefile} --exclude {output.skip_snps} --out {params.baseout} &> {log}
        """

rule create_vcf:
    input: rules.create_plink_lifted.output.plink,
    output: 
        tmp_vcf = temp("output/{project}/liftover/04_hg19_vcf/{project}.{reference}ToHg19.vcf"),
        vcf = "output/{project}/liftover/04_hg19_vcf/{project}.{reference}ToHg19.vcf.gz",
        stats = "output/{project}/liftover/04_hg19_vcf/{project}.{reference}ToHg19.vcf.stats",
    conda: "../envs/liftover.yaml"
    params:
        basein = lambda wildcards,input: path.splitext(input[0])[0],
        baseout =  "output/{project}/liftover/04_hg19_vcf/{project}.{reference}ToHg19",
    log: "logs/liftover/{project}/{project}.{reference}.04_create_vcf.log",
    threads: config['BCFTOOLS_THREADS']
    shell: 
        """
        plink --allow-extra-chr --allow-no-sex --real-ref-alleles --keep-allele-order --snps-only just-acgt --bfile {params.basein} --recode vcf --out {params.baseout} &> {log}
        bcftools +fill-tags {output.tmp_vcf} -- -d  2> {log} | bcftools view -U | bcftools norm --threads {threads} --rm-dup snps -Oz > {output.vcf} 2>> {log}
        bcftools index {output.vcf} &>> {log} 
        bcftools stats {output.vcf} > {output.stats} 2>> {log}
        """
        #bcftools +fill-tags {output.tmp_vcf} -- -d  2> {log} | bcftools view -U | bcftools norm --threads {threads} --rm-dup snps -Oz > {output.vcf} 2>> {log}
        #bcftools +fill-tags {output.tmp_vcf} -- -d  2> {log} | awk '/^#/ {{print}} !($1"_"$2 in a){{a[$1"_"$2]=$1"_"$2;print $0}}' | bcftools view -U -Oz > {output.vcf} 2>> {log}
