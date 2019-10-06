import yaml
import sys
from os import path
from snakemake.utils import validate

# TODO: 
    # 1. Fix shapeit genetic map generator file 
    # 2. Fix shapeit config schema json


configfile: "config.yaml"
validate(config, "config.schema.json")
#report: "report/workflow.rst"

# Allow users to fix the underlying OS via singularity.
singularity: "docker://continuumio/miniconda3"

rule kirimp_ready:
    input:
        [expand('output/{project}/kirimp/02_shapeit/shapeit_v4/{project}.{region}.phased.{out_type}', project=project, out_type = ['haps', 'sample'], region=region)
         for project in config['project'] for region in config['project'][project]['shapeit']['regions']],
        [expand('output/{project}/kirimp/02_shapeit/shapeit_v2/{project}.{region}.phased.{out_type}', project=project, out_type = ['haps', 'sample'], region=region)
         for project in config['project'] for region in config['project'][project]['shapeit']['regions']],

rule kirimp_ready_v4:
    input:
        [expand('output/{project}/kirimp/02_shapeit/shapeit_v2/{project}.{region}.phased.{out_type}', project=project, out_type = ['haps', 'sample'], region=region)
         for project in config['project'] for region in config['project'][project]['shapeit']['regions']],
rule kirimp_panel:
    input:
    	"input/meta/kirimp/%s" % path.basename(config['KIRIMP_PANEL_URL'])

rule liftover:
    input:
    	["input/meta/liftover/%sToHg19.over.chain.gz" % config['project'][project]['liftover']['reference'].lower() for project in config['project']],
        [expand("output/{project}/liftover/04_hg19_vcf/{project}.{reference}ToHg19.vcf.gz", project = project, reference = config['project'][project]['liftover']['reference']) for project in config['project']]

rule kirimp_encode:
    input:
        expand('output/{project}/kirimp/01_freq_encode_snps/{project}.vcf.gz', project=config['project'].keys()),

rule print_defaults:
    run: yaml.dump(config, sys.stdout)
        
rule create_dag:
    output: "dags/dag.svg"
    shell: "snakemake --dag | dot -Tsvg > {output}"
include: "rules/kirimp_panel.smk"
include: "rules/liftover.smk"
include: "rules/freq_encode_snps.smk"
include: "rules/shapeit.smk"
