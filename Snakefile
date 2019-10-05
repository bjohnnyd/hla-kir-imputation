from os import path
from snakemake.utils import validate


configfile: "config.yaml"
validate(config, "config.schema.json")
#report: "report/workflow.rst"

# Allow users to fix the underlying OS via singularity.
singularity: "docker://continuumio/miniconda3"

rule all_kirimp:
    input:
        expand('output/{project}/kirimp/01_freq_encode_snps/{project}.vcf.gz', project=config['project'].keys())

rule kirimp_panel:
    input:
    	"input/meta/kirimp/%s" % path.basename(config['KIRIMP_PANEL_URL'])

rule liftover:
    input:
    	["input/meta/liftover/%sToHg19.over.chain.gz" % config['project'][project]['liftover']['reference'].lower() for project in config['project']],
        [expand("output/{project}/liftover/04_hg19_vcf/{project}.{reference}ToHg19.vcf.gz", project = project, reference = config['project'][project]['liftover']['reference']) for project in config['project']]

include: "rules/kirimp_panel.smk"
include: "rules/liftover.smk"
include: "rules/freq_encode_snps.smk"
