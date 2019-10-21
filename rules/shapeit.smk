import re
from os import path
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider


HTTP = HTTPRemoteProvider()

rule generate_kirimp_gmap:
    input: HTTP.remote(config['SHAPEIT_GENMAP_URL'], keep_local=False, allow_redirects=True)
    output: 'input/meta/shapeit/kirimp.chr19.gmap.txt.gz'
    log : 'logs/generate_kirimp_gmap/kirimp.chr19.gmap.log'
    shell: "cat {input} | tar xzfO - --wildcards 'chr19.b37.gmap.gz' 2> {log}  > {output}"

rule shapeit4:
    input: 
        vcf = lambda wc: config['project'][wc.project]['shapeit']['vcf'],
        gmap = lambda wc: config['project'][wc.project]['shapeit']['gmap'],
    output:
        phased_vcf = 'output/{project}/kirimp/02_shapeit/shapeit_v4/{project}.{region}.phased.vcf.gz',
        filtered_vcf = 'output/{project}/kirimp/02_shapeit/shapeit_v2/{project}.{region}.filtered.vcf.gz',
        hap = 'output/{project}/kirimp/02_shapeit/shapeit_v4/{project}.{region}.phased.hap',
        sample = 'output/{project}/kirimp/02_shapeit/shapeit_v4/{project}.{region}.phased.samples',
        legend = temp('output/{project}/kirimp/02_shapeit/shapeit_v4/{project}.{region}.phased.legend.gz'),
        log = 'output/{project}/kirimp/02_shapeit/shapeit_v4/{project}.{region}.shapeit.log',
    conda: "../envs/shapeit.yml"
    log:
    	shapeit4= "logs/shapeit/{project}/{project}.{region}.shapeit4.err",
    	bcftools= "logs/shapeit/{project}/{project}.{region}.bcftools.shapeit4.err",
    params: 
    	pbwt=lambda wc: config['project'][wc.project]['shapeit']['pbwt'],
    	pbwt_modulo=lambda wc: config['project'][wc.project]['shapeit']['pbwt-modulo'],
    	additional=lambda wc: config['project'][wc.project]['shapeit']['v4_additional'],
        missing_threshold=lambda wc: config['project'][wc.project]['shapeit']['min_missing'],
        out='output/{project}/kirimp/02_shapeit/shapeit_v4/{project}.{region}.phased'
    threads: config['SHAPEIT_THREADS']
    shell: 
        """
            bcftools filter -r {wildcards.region} -Oz -e 'MISS > {params.missing_threshold}' {input.vcf} > {output.filtered_vcf} 2> {log.bcftools}
            bcftools index -f {output.filtered_vcf} 2>> {log.bcftools}
            shapeit4 --thread {threads} --input {output.filtered_vcf} --map {input.gmap} --output {output.phased_vcf} \
            --region {wildcards.region} --pbwt-depth {params.pbwt} --pbwt-modulo {params.pbwt_modulo} --log {output.log} {params.additional} 2> {log.shapeit4}
            bcftools convert --hapsample {params.out} {output.phased_vcf} 2>> {log.bcftools}
            gunzip {output.hap}.gz 
        """

rule shapeit:
    input: 
        vcf = lambda wc: config['project'][wc.project]['shapeit']['vcf'],
        gmap = lambda wc: config['project'][wc.project]['shapeit']['gmap']
    output:
        filtered_vcf = 'output/{project}/kirimp/02_shapeit/shapeit_v2/{project}.{region}.filtered.vcf.gz',
        haps = 'output/{project}/kirimp/02_shapeit/shapeit_v2/{project}.{region}.phased.haps',
        sample = 'output/{project}/kirimp/02_shapeit/shapeit_v2/{project}.{region}.phased.sample',
        graph = 'output/{project}/kirimp/02_shapeit/shapeit_v2/{project}.{region}.graph',
        log = 'output/{project}/kirimp/02_shapeit/shapeit_v2/{project}.{region}.shapeit.log',
    conda: "../envs/shapeit.yml"
    log:
    	shapeit= "logs/shapeit/{project}/{project}.{region}.shapeit.err",
    	bcftools= "logs/shapeit/{project}/{project}.{region}.bcftools.shapeit.err",
    params: 
    	states=lambda wc: config['project'][wc.project]['shapeit']['states'],
        missing_threshold=lambda wc: config['project'][wc.project]['shapeit']['min_missing'],
    	additional=lambda wc: config['project'][wc.project]['shapeit']['v2_additional'],
    threads: config['SHAPEIT_THREADS']
    shell: 
        "bcftools filter -r {wildcards.region} -Oz -e 'MISS > {params.missing_threshold}' {input.vcf} > {output.filtered_vcf} 2> {log.bcftools} && "
        "bcftools index -f {output.filtered_vcf} 2> {log.bcftools}  && "
        "shapeit --thread {threads} --input-vcf {output.filtered_vcf} -M {input.gmap} --states {params.states} "
        "-O {output.haps} {output.sample} --output-graph {output.graph} "
        "--output-log {output.log} {params.additional} 2> {log.shapeit}"
