import re
from os import path
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider

FTP = FTPRemoteProvider()

rule generate_kirimp_gmap:
    input: FTP.remote(config['KIRIMP_GENMAP_URL'], keep_local=False)
    output: 'input/meta/shapeit/kirimp.chr19.gmap.txt.gz'
    log : 'logs/generate_kirimp_gmap/kirimp.chr19.gmap.log'
    shell: "cat {input} | tar xzfO - --wildcards '*_chr19.txt' 2> {log} | sed -E 's:^chr::' | gzip > {output}"

rule shapeit4:
    input: 
        vcf = lambda wc: config['project'][wc.project]['shapeit']['vcf'],
        gmap = lambda wc: config['project'][wc.project]['shapeit']['gmap'],
    output:
        phased_vcf = 'output/{project}/kirimp/02_shapeit/shapeit_v4/{project}.{region}.phased.vcf.gz',
        haps = 'output/{project}/kirimp/02_shapeit/shapeit_v4/{project}.{region}.phased.haps',
        sample = 'output/{project}/kirimp/02_shapeit/shapeit_v4/{project}.{region}.phased.sample',
        graph = 'output/{project}/kirimp/02_shapeit/shapeit_v4/{project}.{region}.graph',
        log = 'output/{project}/kirimp/02_shapeit/shapeit_v4/{project}.{region}.shapeit.log',
    conda: "../envs/shapeit.yml"
    log: "logs/shapeit/{project}/{project}.{region}.shapeit4.err"
    params: 
    	pbwt=lambda wc: config['project'][wc.project]['shapeit']['pbwt'],
    	pbwt_modulo=lambda wc: config['project'][wc.project]['shapeit']['pbwt-modulo'],
    	additional=lambda wc: config['project'][wc.project]['shapeit']['v4_additional'],
        missing_threshold=lambda wc: config['project'][wc.project]['shapeit']['min_missing'],
    threads: config['SHAPEIT_THREADS']
    shell: 
        """
            bcftools filter -r {wildcards.region} -e 'MISS > {params.missing_threshold}' {input.vcf} | 
            shapeit4 --thread {threads} --input - --map {input.gmap} --output {output.phased_vcf} \
            --region {wildcards.region} --pbwt-depth {params.pbwt} --pbwt-modulo {params.pbwt_modulo} --log {output.log} {params.additional} 2> {log}
            bcftools convert -hapsample {output.haps},{output.sample} {output.phased_vcf}
        """

rule shapeit:
    input: 
        vcf = lambda wc: config['project'][wc.project]['shapeit']['vcf'],
        gmap = lambda wc: config['project'][wc.project]['shapeit']['gmap']
    output:
        haps = 'output/{project}/kirimp/02_shapeit/shapeit_v2/{project}.{region}.phased.haps',
        sample = 'output/{project}/kirimp/02_shapeit/shapeit_v2/{project}.{region}.phased.sample',
        graph = 'output/{project}/kirimp/02_shapeit/shapeit_v2/{project}.{region}.graph',
        log = 'output/{project}/kirimp/02_shapeit/shapeit_v2/{project}.{region}.shapeit.log',
    conda: "../envs/shapeit.yml"
    log: "logs/shapeit/{project}/{project}.{region}.shapeit.err"
    params: 
    	states=lambda wc: config['project'][wc.project]['shapeit']['states'],
        missing_threshold=lambda wc: config['project'][wc.project]['shapeit']['min_missing'],
    	additional=lambda wc: config['project'][wc.project]['shapeit']['v2_additional'],
    threads: config['SHAPEIT_THREADS']
    shell: 
        "bcftools filter -r {wildcards.region} -e 'MISS > {params.missing_threshold}' {input.vcf} | "
        "shapeit --threads {threads} --input-vcf - -M {input.gmap} --states {params.states} "
        "-O {output.haps} {output.sample} --output-graph {output.graph} "
        "--output-log {output.log} {params.additional} 2> {log}"
