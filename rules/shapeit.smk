import re
from os import path

from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()

# curl -L ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/working/20110106_recombination_hotspots/HapmapII_GRCh37_RecombinationHotspots.tar.gz --output - | tar xzfv - --wildcards '*_chr19.txt'

def get_gmap(wildcards):
    fpath = config['project'][wildcards.project]['shapeit']['gmap']
    if re.match("^(ftp|http)",fpath):
        HTTP.remote(fpath, keep_local=True)

    

rule generate_gmap:
    input: get_gmap
    output: 'output/{project}/02_gmap/{project}.gmap'
    params: '--wildcards %s' ' '.join["'*_chr%d.txt'" for chromosome in config['project'][wc.project]['shapeit']['chromosomes']  
    shell:
        if remote and ends in tar:
            shell(
                "cat {input} --output - | tar xOzf {params} | sed -E 's:^chr::g' > {output}"
            )
        else:
            shell(
                "mv {input} {output}"
            )

rule shapeit:
    input: 
        vcf = rules.freq_encode_snps.output.vcf,
        gmap = rules.generate_gmap.output[0]
    output:
        haps = 'output/{project}/kirimp/02_shapeit/{project}.phased.haps.gz',
        samples = 'output/{project}/kirimp/02_shapeit/{project}.phased.sample.gz',
        graph = 'output/{project}/kirimp/02_shapeit/{project}.graph',
        log = 'output/{project}/kirimp/02_shapeit/{project}.shapeit.log',
    conda: "../envs/shapeit.yml"
    log: "logs/shapeit/{project}/{project}.err"
    params: lambda wc: config['project'][wc.project]['shapeit']['additional']
    shell: 
        "shapeit --input-vcf {input.vcf} -M {input.gmap} "
        "-O {output.haps} {output.sample} --output-graph {output.graph} "
        "--output-log {output.log} {additional} 2> {log}"