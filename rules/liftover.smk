from os import path
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()

URL_DICT = {
    "HG18" : "http://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz",
    "HG17" : "https://hgdownload.soe.ucsc.edu/goldenPath/hg17/liftOver/hg17ToHg19.over.chain.gz",
    "HG16" : "https://hgdownload.soe.ucsc.edu/goldenPath/hg16/liftOver/hg16ToHg19.over.chain.gz",
}

def get_overchain_url(wildcards):
    reference=config["project"][wildcards.project]['liftover']['reference']
    overchain_url=config['OVERCHAIN_URLS'][reference]
    print(HTTP.remote(overchain_url, keep_local=True))
    HTTP.remote(overchain_url, keep_local=True)


rule get_overchain:
    input: get_overchain_url   
    output: path.join("input/meta/liftover/",path.basename("over/chain.txt"))
    shell: "mv {input} {output}"

