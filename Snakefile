from os import path
from snakemake.utils import validate


configfile: "config.yaml"
validate(config, "config.schema.yaml")
#report: "report/workflow.rst"

# Allow users to fix the underlying OS via singularity.
singularity: "docker://continuumio/miniconda3"

rule all_kirimp:
    input:
    	"input/meta/%s" % path.basename(config['KIRIMP_PANEL_URL'])

rule all_bed_to_b37:
    input:
    	"input/meta/%s" % path.basename(config['KIRIMP_PANEL_URL'])

include: "rules/kirimp_panel.smk"
