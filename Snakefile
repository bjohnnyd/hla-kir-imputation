from os import path
from snakemake.utils import validate
# The main entry point of your workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.


configfile: "config.yaml"
validate(config, "config.schema.yaml")
#validate(config, "config.schema.json")
#report: "report/workflow.rst"

# Allow users to fix the underlying OS via singularity.
#singularity: "docker://continuumio/miniconda3"

rule all_kirimp:
    input:
    	"input/meta/%s" % path.basename(config['KIRIMP_PANEL_URL'])


include: "rules/kirimp_panel.smk"
