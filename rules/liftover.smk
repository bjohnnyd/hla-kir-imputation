from os import path
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
from snakemake.utils import format as snakeformat

HTTP = HTTPRemoteProvider()

URL_DICT = {
    "HG18" : "http://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz",
    "HG17" : "https://hgdownload.soe.ucsc.edu/goldenPath/hg17/liftOver/hg17ToHg19.over.chain.gz",
    "HG16" : "https://hgdownload.soe.ucsc.edu/goldenPath/hg16/liftOver/hg16ToHg19.over.chain.gz",
}

rule get_overchain:
    input: lambda wc: HTTP.remote(URL_DICT[wc.reference.upper()], keep_local = True)
    output: "input/meta/liftover/{reference}ToHg19.over.chain.gz"
    shell: "mv {input} {output}"

rule plink_to_ucsc:
    input:
        plink = lambda wc: config['project'][wc.project]['liftover']['plink'],
        chain = "input/meta/liftover/{reference}ToHg19.over.chain.gz",
    output:
        plink_map = "output/{project}/liftover/{reference}/{project}.{reference}.map",
        ucsc_bed = "output/{project}/liftover/{reference}/{project}.{reference}.ucsc.bed",
        plink_changefile = "output/{project}/liftover/hg19/{project}.{reference}.changefile",
        plink_changefile_unlifted = "output/{project}/liftover/hg19/{project}.{reference}.changefile.unlifted",
        hg19_unlifted = "output/{project}/liftover/hg19/{project}.{reference}ToHg19.unlifted",
        hg19_plink = "output/{project}/liftover/hg19/{project}.{reference}ToHg19.plink.bed",
        hg19_ucsc_bed = "output/{project}/liftover/hg19/{project}.{reference}ToHg19.ucsc.bed",
        hg19_vcf = "output/{project}/liftover/hg19/{project}.{reference}ToHg19.vcf.gz",
    params:
        plink_in = lambda widcards, output, input: path.splitext(input.plink)[0],
        plink_out = "output/{project}/liftover/{reference}/{project}.{reference}",
        plink_lifted_out = "output/{project}/liftover/hg19/{project}.{reference}ToHg19",
    conda: "../envs/liftover.yaml"
    log: 
    	plink = "logs/liftover/{project}/{project}.{reference}.plink.err",
        liftover =  "logs/liftover/{project}/{project}.{reference}.liftOver.err",
        bcftools =  "logs/liftover/{project}/{project}.{reference}.bcftools.err"
    shell: 
        r"""
        plink --bfile {params.plink_in} --allow-no-sex --real-ref-alleles --recode --out {params.plink_out} &> {log.plink} && 
        awk 'BEGIN{{OFS="\t"}} {{print "chr"$1,$4-1,$4,$2,$3}}' {output.plink_map} > {output.ucsc_bed} &&
        liftOver {output.ucsc_bed} {input.chain} {output.hg19_ucsc_bed} {output.hg19_unlifted} &> {log.liftover} &&
        awk '{{print $2,$4}}' {output.hg19_ucsc_bed} > {output.plink_changefile} &&
        awk '{{print $2,$4}}' {output.hg19_unlifted} > {output.plink_changefile_unlifted} &&
        plink --file {params.plink_out} --allow-no-sex --real-ref-alleles --make-bed --update-map {output.plink_changefile} --exclude {output.plink_changefile_unlifted} --out {params.plink_lifted_out}.plink &>> {log.plink} &&
        plink --allow-extra-chr --allow-no-sex --real-ref-alleles --keep-allele-order --snps-only just-acgt --bfile {params.plink_lifted_out}.plink --recode vcf --out {params.plink_lifted_out} &>> {log.plink} &&
        bcftools +fill-tags {params.plink_lifted_out}.vcf -Oz -- -d  2> {log.bcftools} |   bcftools annotate -Oz --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' > {output.hg19_vcf} 2>> {log.bcftools}  && bcftools index {output.hg19_vcf} &>> {log.bcftools} 
        """
