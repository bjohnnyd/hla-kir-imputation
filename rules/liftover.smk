from os import path
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

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