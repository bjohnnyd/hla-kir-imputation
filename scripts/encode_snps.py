#!/usr/bin/env python
import argparse
from os import path
from cyvcf2 import VCF, Writer


CHROMOSOME_19_ANNOTATION = {
    'ucsc': 'chr19',
    'ensembl': '19',
    'genbank': 'CM000681.2'
}


def main(arguments=None):
    args = parse_arguments()
    vcf = VCF(args['vcf_file'])
    panel = generate_panel_data(panel_file=args['reference_panel'],
                                chr=args['chromosomes'], annotation=args['chromosome_annotation'],
                                panel_type=args['reference_panel_type'], panel_format=args['reference_panel_format'])
    print(panel)
    print(vcf.seqnames)


def generate_panel_data(panel_file, chr=None, annotation='ensembl', panel_type='kirimp', panel_format='id,pos,allele0,allele1,allele1_frequency'):
    if panel_type == 'kirimp':
        snp_dict = kirimp_parser(
            kirimp=panel_file, chr=CHROMOSOME_19_ANNOTATION[annotation])
        return(snp_dict)
    if panel_type == 'custom':


def kirimp_parser(kirimp, chr):
    snp_dict = {}
    with open(kirimp) as f:
        next(f)
        for line in f:
            snp = line.split(',')
            snp_dict[chr+'_'+snp[1]] = {
                "A1": snp[2],
                "A2": snp[3],
                "freq": float(snp[4].strip())
            }
    return(snp_dict)


def parse_arguments(arguments=None):
    parser = argparse.ArgumentParser(
        description="This script encodes SNPs in a VCF to a reference panel based on allele frequencies",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        "-v", "--vcf-file", help="VCF/BCF file to re-encode (can be compressed with bgzip)", required=True, type=str)
    parser.add_argument(
        "-r", "--reference-panel", help="Reference panel file  containing data in the format [chrom pos freq] or [pos freq] if only one chromosome in input VCF", required=True, type=str)
    parser.add_argument(
        "-rt", "--reference-panel-type", help="Reference panel file type", choices=["kirimp", "custom"], default='kirimp', type=str)
    parser.add_argument(
        "-rf", "--reference-panel-format", help="Custom reference panel format type", type=str)
    parser.add_argument(
        "-o", "--output", help="Output vcf file", required=False, default='stdout', type=str)
    parser.add_argument(
        "-O", "--output-type", help="Output vcf file type", choices=["z", "v", "b"], type=str, default='z')
    parser.add_argument(
        "-chr", "--chromosomes", help="Chromosome over which to encode SNPs ", required=False, nargs='?', type=str)
    parser.add_argument(
        "--chromosome-annotation", help="Chromosome annotation type in the VCF", choices=['ucsc', 'ensembl', 'genbank'], default='ensembl', type=str)
    parser.add_argument("-a", "--ambigious",
                        help="Determines whether ambigious alternate alleles should be dropped", action='store_false')
    parser.add_argument("-c", "--fix-complement-ref-alt",
                        help="Should ref/alt that are complements be fixed with respect to frequency", action='store_false')
    parser.add_argument("-min", "--min-ambigious-threshold", help="Alternate alleles above this frequency and below the max ambigious frequency will be flagged as ambigious",
                        default=0.49, type=float)
    parser.add_argument("-max", "--max-ambigious-threshold", help="Alternate alleles above this frequency and below the max ambigious frequency will be flagged as ambigious",
                        default=0.51, type=float)
    args = vars(parser.parse_args())
    print(args)
    if args['reference_panel_type'] == 'custom' and args['reference_panel_format'] is None:
        parser.error(
            'custom --reference-panel-type requires --reference-panel-format to be set')
    return(args)


if __name__ == '__main__':
    main()
