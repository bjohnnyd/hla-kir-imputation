#!/usr/bin/env python
import argparse
from os import path
from cyvcf2 import VCF, Writer


def main(arguments=None):
    args = parse_arguments()
    vcf = VCF(args['vcf-file'])
    print(vcf.seqnames)


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
    parser.add_argument("-a", "--ambigious",
                        help="Determines whether ambigious alternate alleles should be dropped", action='store_false')
    parser.add_argument("-c", "--fix-complement-ref-alt",
                        help="Should ref/alt that are complements be fixed with respect to frequency", action='store_false')
    parser.add_argument("-min", "--min-ambigious-threshold", help="Alternate alleles above this frequency and below the max ambigious frequency will be flagged as ambigious",
                        default=0.49, type=float)
    parser.add_argument("-max", "--max-ambigious-threshold", help="Alternate alleles above this frequency and below the max ambigious frequency will be flagged as ambigious",
                        default=0.51, type=float)
    args = vars(parser.parse_args())
    if args['reference-panel-type'] == 'custom' and args['reference-panel-format'] is None:
        parser.error(
            'custom --reference-panel-type requires --reference-panel-format to be set')
    return(args)

if __name__ == '__main__':
    main()
