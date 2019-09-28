from os import path
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()

rule get_kirimp_panel:
    input: HTTP.remote(config["KIRIMP_PANEL_URL"], keep_local=True)
    output: path.join("input/meta/kirimp",path.basename(config["KIRIMP_PANEL_URL"]))
    shell: "mv {input} {output}"

