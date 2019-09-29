#!/usr/bin/env python
import argparse
from os import path
from cyvcf2 import VCF, Writer


CHROMOSOME_19_ANNOTATION = {
    'ucsc': 'chr19',
    'ensembl': '19',
    'genbank': 'CM000681.2'
}

KIRIMP_HEADER = [
    "id",
    "position",
    "allele0",
    "allele1",
    "allele1_frequency"
]

CUSTOM_HEADER = [
    "chrom",
    "pos",
    "a0",
    "a1",
    "freq"
]


""" Class to represent genotypes in 0/1 format might not be necessary as I can flip from there"""


class Genotype(object):
    __slots__ = ('alleles', 'phased')

    def __init__(self, li):
        self.alleles = li[:-1]
        self.phased = li[-1]

    def __str__(self):
        sep = "/|"[int(self.phased)]
        return sep.join("0123."[a] for a in self.alleles)
    __repr__ = __str__


"""
Custom panel format needs to be: chrom,pos,a0,a1,freq
"""


def main(arguments=None):
    args = parse_arguments()
    vcf = VCF(args['vcf_file'])
    panel = generate_panel_data(
        panel_file=args['reference_panel'], chr=args['chromosomes'], annotation=args['chromosome_annotation'], panel_type=args['reference_panel_type'])
    for variant in vcf:
        variant_id_start = str(variant.CHROM) + '_' + str(variant.start)
        variant_id_end = str(variant.CHROM) + '_' + str(variant.end)
        if variant_id_start in panel or variant_id_end in panel:
            print(variant.INFO.get('AF'))
            print(variant.genotypes)
            print(variant.gt_types)


def generate_panel_data(panel_file, chr=None, annotation='ensembl', panel_type='kirimp'):
    f = open(panel_file, 'r')
    header_line = next(f).strip()
    sep = get_separator(header_line)
    header_line = [cell.replace('"', '')
                   for cell in header_line.split(sep)]
    if panel_type == 'kirimp':
        chromosome = CHROMOSOME_19_ANNOTATION[annotation]
        if header_line != KIRIMP_HEADER:
            raise TypeError(
                "If input panel type is kirimp, the panel needs to contain a comma-separated header:\n%s" % ','.join(KIRIMP_HEADER))
    else:
        if header_line != CUSTOM_HEADER:
            raise TypeError(
                "If input panel type is custom, the panel needs to contain a comma-separated header:\n%s" % ','.join(CUSTOM_HEADER))
    snp_dict = {
        (chromosome+'_'+cells[1] if panel_type == 'kirimp' else cells[0]+'_'+cells[1]): {
            "A0": cells[2],
            "A1": cells[3],
            "freq": float(cells[4].strip())
        } if panel_type == 'kirimp'
        else {
            'A0': cells[2],
            'A1': cells[3],
            'freq': float(cells[4])
        } for cells in [line.strip().replace('"', '').split(sep) for line in f]}
    f.close()
    return(snp_dict)


def get_separator(line):
    tabs = line.count(r'\t')
    commas = line.count(r',')
    if tabs == 4 and commas != 4:
        sep = r'\t'
    elif tabs != 4 and commas == 4:
        sep = ','
    else:
        raise TypeError(
            'Cannot determine separator from file please specify separator directly as an argument [--reference-panel-col-separator]')
    return(sep)


# def kirimp_parser(kirimp, chr, type):
#     snp_dict = {}
#     with open(kirimp) as f:
#         header = [cell.replace(''"','') for cell in next(f).strip().split(', ')]
#         print(header)
#         if header != KIRIMP_HEADER:
#             raise TypeError(
#                 "If input panel type is kirimp, the panel needs to contain a comma-separated header:\n%s" % ','.join(KIRIMP_HEADER))
#         for line in f:
#             snp = line.split(',')
#             snp_dict[chr+'_'+snp[1]] = {
#                 "A1": snp[2],
#                 "A2": snp[3],
#                 "freq": float(snp[4].strip())
#             }
#     return(snp_dict)


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
    if args['reference_panel_type'] == 'custom' and args['reference_panel_format'] is None:
        parser.error(
            'custom --reference-panel-type requires --reference-panel-format to be set')
    return(args)


if __name__ == '__main__':
    main()
