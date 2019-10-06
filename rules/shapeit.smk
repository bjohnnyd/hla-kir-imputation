import re
from os import path
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider

FTP = FTPRemoteProvider()

rule generate_kirimp_gmap:
    input: FTP.remote(config['KIRIMP_GENMAP_URL'], keep_local=False)
    output: 'output/kirimp/gmap/chr19/kirimp.chr19.gmap.txt.gz'
    log : 'logs/generate_kirimp_gmap/kirimp.chr19.gmap.log'
    shell: "cat {input} | tar xzfO - --wildcards '*_chr19.txt' 2> {log} | sed -E 's:^chr::' | gzip > {output}"

rule shapeit:
    input: 
        vcf = rules.freq_encode_snps.output.vcf,
        gmap = lambda wc: config['project'][wc.project]['shapeit']['gmap']
    output:
        phased_vcf = 'output/{project}/kirimp/02_shapeit/{project}.phased.vcf.gz',
        haps = 'output/{project}/kirimp/02_shapeit/{project}.phased.haps.gz',
        sample = 'output/{project}/kirimp/02_shapeit/{project}.phased.sample.gz',
        graph = 'output/{project}/kirimp/02_shapeit/{project}.graph',
        log = 'output/{project}/kirimp/02_shapeit/{project}.shapeit.log',
    conda: "../envs/shapeit.yml"
    log: "logs/shapeit/{project}/{project}.err"
    params: 
    	pbwt=lambda wc: config['project'][wc.project]['shapeit']['pbwt'],
    	pbwt-modulo=lambda wc: config['project'][wc.project]['shapeit']['pbwt-modulo'],
    	region=lambda wc: config['project'][wc.project]['shapeit']['regions'],
    	additonal=lambda wc: config['project'][wc.project]['shapeit']['additional'],
    threads: config['SHAPEIT_THREADS']
    shell: 
	"""
        shapeit4 --thread {threads} --input {input.vcf} --map {input.gmap} --output {output.phased_vcf} \
        --region {params.region} --pbwt-depth {params.pbwt} --log {output.log} {params.additional} 2> {log}
	"""
#        "shapeit --threads {threads} --input-vcf {input.vcf} -M {input.gmap} "
#        "-O {output.haps} {output.sample} --output-graph {output.graph} "
#        "--output-log {output.log} {params} 2> {log}"
