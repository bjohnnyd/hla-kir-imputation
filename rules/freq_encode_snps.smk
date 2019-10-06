from os import path

def generate_vcf_input(wildcards):
    input_vcf = config['project'][wildcards.project]['freq_encode_snps']['vcf']
    if input_vcf == rules.create_vcf.output.vcf and '{reference}' in input_vcf:
        reference = config['project'][wildcards.project]['liftover']['reference']
        return input_vcf.replace('{reference}','%s') % reference.lower()
    else:
        return input_vcf

rule freq_encode_snps:
    input: 
        panel = "input/meta/kirimp/%s" % path.basename(config['KIRIMP_PANEL_URL']),
        vcf = generate_vcf_input,
        # vcf = lambda wc: config['project'][wc.project]['freq_encode_snps']['vcf'],
    output:
        vcf = 'output/{project}/kirimp/01_freq_encode_snps/{project}.vcf.gz',
        index = 'output/{project}/kirimp/01_freq_encode_snps/{project}.vcf.gz.csi',
        plot = report('output/{project}/kirimp/01_freq_encode_snps/{project}.png', caption = 'report/freq_encode_snps.rst', category = "Frequency Encode SNPs to KIR*IMP Panel"),
    conda: "../envs/freq_encode_snps.yml"
    log: "output/{project}/kirimp/01_freq_encode_snps/{project}_log.toml"
    params: lambda wc: config['project'][wc.project]['freq_encode_snps']['additional']
    shell: 
        "scripts/frequency_encode_snps.py -v {input.vcf} "
        "-r {input.panel} -o {output.vcf} {params} 2> {log} && "
        "bcftools index {output.vcf}"