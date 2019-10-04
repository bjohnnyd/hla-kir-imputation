#!/usr/bin/env python
import argparse
import numpy as np
from os import path
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
from cyvcf2 import VCF, Writer


plt.rcParams["figure.figsize"] = (16, 12)

CHROMOSOME_19_ANNOTATION = {"ucsc": "chr19", "ensembl": "19", "genbank": "CM000681.2"}

# Required headers for panel input files, either kirimp (kirimp.uk1.snp.info.csv) or custom with the required fields
KIRIMP_HEADER = ["id", "position", "allele0", "allele1", "allele1_frequency"]
CUSTOM_HEADER = ["chrom", "pos", "a0", "a1", "freq"]

COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G", ".": "."}

""" CUSTOM HEADERS """
UPDATED = {
    "ID": "UPD",
    "Description": "REF and ALT updated based on reference panel frequency",
    "Type": "Flag",
    "Number": "A",
}

PANEL_FREQ_DIF = {
    "ID": "PFD",
    "Description": "Alternate frequency difference to reference panel frequency",
    "Type": "Float",
    "Number": "A",
}

PANEL_FREQ_DIFF = {
    "ID": "PFD",
    "Description": "Alternate frequency difference to reference panel frequency",
    "Type": "Float",
    "Number": "A",
}

MISSINGNES = {
    "ID": "MISS",
    "Description": "Missing Genotype Frequency",
    "Type": "Float",
    "Number": "A",
}

MAF = {
    "ID": "MAF",
    "Description": "Minor Allele Frequency",
    "Type": "Float",
    "Number": "A",
}

""" Class to represent genotypes in 0/1 format might not be necessary as I can flip from there"""


class Genotype(object):
    __slots__ = ("alleles", "phased")
    FLIP_DICT = {0: 1, 1: 0, -1: -1}

    def __init__(self, li):
        self.alleles = li[:-1]
        self.phased = li[-1]

    def __str__(self):
        sep = "/|"[int(self.phased)]
        return sep.join("0123."[a] for a in self.alleles)

    def flip(self):
        self.alleles = [Genotype.FLIP_DICT[allele] for allele in self.alleles]

    def genotype(self):
        return self.alleles + [self.phased]

    __repr__ = __str__


""" Class to keep track of VCF file summaries for variants """


class VCFSummary(object):
    __slots__ = ("ambigious", "unknown_alt", "updated", "flipped", "VARIANTS")

    def __init__(self):
        self.VARIANTS = {}
        self.ambigious = 0
        self.unknown_alt = 0
        self.updated = 0
        self.flipped = 0

    def add_ambigious(self):
        self.ambigious += 1

    def add_unknown_alt(self):
        self.unknown_alt += 1

    def add_updated(self):
        self.updated += 1

    def add_flipped(self):
        self.flipped += 1

    def add_variant(self, v_id):
        self.VARIANTS[v_id] = {"freq": None, "updated_freq": None, "panel_freq": None}

    def add_variant_dict(self, vdict):
        self.VARIANTS.update(vdict)


def main(arguments=None):
    args = parse_arguments()
    vcf = VCF(args["vcf_file"])
    vcf.add_info_to_header(UPDATED)
    vcf.add_info_to_header(PANEL_FREQ_DIFF)
    vcf.add_info_to_header(MISSINGNES)
    vcf.add_info_to_header(MAF)

    w = Writer(args["output"], vcf)
    panel = generate_panel_data(
        panel_file=args["reference_panel"],
        chr=args["chromosomes"],
        annotation=args["chromosome_annotation"],
        panel_type=args["reference_panel_type"],
        separator=args["separator"],
    )

    vcf_summary = VCFSummary()

    for variant in vcf:
        variant_id_end = str(variant.CHROM) + "_" + str(variant.end)
        if variant_id_end in panel:
            panel_variant = panel[variant_id_end]
            if not variant.ALT:
                print("-" * 100)
                print("No alternate called or missing for variant %s" % variant.ID)
                vcf_summary.add_unknown_alt()
                continue
            if (
                args["ambigious"]
                and variant.aaf > args["min_ambigious_threshold"]
                and variant.aaf < args["max_ambigious_threshold"]
            ):
                vcf_summary.add_ambigious()
                continue
            if should_recode(variant, panel_variant):
                swap_ref_alt(variant)
                variant.INFO["UPD"] = True
                vcf_summary.add_updated()
            if (
                should_flipstrand(variant, panel_variant)
                and args["fix_complement_ref_alt"]
            ):
                flipstrand(variant)
                variant.INFO["UPD"] = True
                vcf_summary.add_flipped()

            vcf_summary.add_variant(variant_id_end)
            v_freq = variant.INFO.get("AF")

            variant.INFO["PFD"] = abs(variant.INFO.get("AF") - panel_variant["freq"])
            variant.INFO["MISS"] = np.sum(variant.gt_types == 2) / len(variant.gt_types)
            variant.INFO["MAF"] = v_freq if v_freq < 0.5 else 1 - v_freq

            vcf_summary.VARIANTS[variant_id_end].update(
                {
                    "freq": variant.aaf,
                    "panel_freq": panel_variant["freq"],
                    "updated_freq": v_freq,
                    "MAF": variant.INFO.get("MAF"),
                    "MISS": variant.INFO.get("MISS"),
                    "PFD": variant.INFO.get("PFD"),
                }
            )
        w.write_record(variant)
    w.close()
    vcf.close()
    create_summary_plot(vcf_summary, outfile=args["output"].split(".")[0] + ".png")


def swap_ref_alt(variant):
    gts = variant.genotypes
    gts = [Genotype(li) for li in gts]
    for gt in gts:
        gt.flip()
    variant.genotypes = [gt.genotype() for gt in gts]
    updated_frequency = sum([gt.alleles.count(1) for gt in gts]) / (
        2 * len([gt for gt in gts if -1 not in gt.alleles])
    )
    temp_nuc = variant.REF
    variant.REF = variant.ALT[0]
    variant.ALT = [temp_nuc]
    variant.INFO["AF"] = updated_frequency


def flipstrand(variant, COMPLEMENT=COMPLEMENT):
    variant.REF = COMPLEMENT[variant.REF]
    variant.ALT = COMPLEMENT[variant.ALT]
    swap_ref_alt(variant)


def should_recode(variant, panel_variant):
    panel_nucleotides = [panel_variant["A0"], panel_variant["A1"]]
    variant_nucleotides = variant.ALT[:]
    variant_nucleotides.extend(variant.REF)
    frequency_synced = (
        panel_variant["freq"] > 0.5 and variant.INFO.get("AF") > 0.5
    ) or (panel_variant["freq"] < 0.5 and variant.INFO.get("AF") < 0.5)
    nucleotides_synced = all(nuc in variant_nucleotides for nuc in panel_nucleotides)
    return not (frequency_synced and nucleotides_synced)


def should_flipstrand(variant, panel_variant, COMPLEMENT=COMPLEMENT):
    unsynced = should_recode(variant, panel_variant)
    is_alt_complement = COMPLEMENT[variant.REF] == variant.ALT
    return unsynced and is_alt_complement


def create_summary_plot(v_summary, outfile):
    X_OFFSET = 0.01
    Y_OFFSET = 0.2

    v_ids = [v_id for v_id in v_summary.VARIANTS]

    fig = plt.figure()

    ax1 = fig.add_subplot(321)
    ax2 = fig.add_subplot(323, sharex=ax1)
    ax3 = fig.add_subplot(325)
    ax4 = fig.add_subplot(222)
    ax5 = fig.add_subplot(224, sharex=ax4)

    (freq, panel_freq, updated_freq) = zip(
        *[
            (
                v_summary.VARIANTS[v_id]["freq"],
                v_summary.VARIANTS[v_id]["panel_freq"],
                v_summary.VARIANTS[v_id]["updated_freq"],
            )
            for v_id in v_ids
        ]
    )

    coefs = (
        np.corrcoef(freq, panel_freq)[1, 0],
        np.corrcoef(updated_freq, panel_freq)[1, 0],
    )
    titles = (
        "Original VCF Frequencies Compared to Panel Frequencies",
        "Updated VCF Frequencies Compared to Panel Frequencies",
    )

    freqs = (freq, updated_freq)

    for i, ax in enumerate([ax1, ax2]):
        ax.set_title(titles[i], fontsize=9)
        ax.scatter(freqs[i], panel_freq, s=10, alpha=0.7)
        ax.annotate(
            "corr = %.2f" % coefs[i],
            (max(freqs[i]) - X_OFFSET, max(panel_freq) - Y_OFFSET),
            ha="center",
        )
        ax.set_ylabel("Panel Allele Frequency", fontsize=8)
        ax.set_xlabel("VCF Alternate Frequency", fontsize=8)
        ax.plot(ax.get_xlim(), ax.get_ylim(), ls="--", c=".3")
        # for i, vid in enumerate(ids):
        #     if abs(updated_freq[i] - panel_freq[i]) > 0.1:
        #         print(vid)
        #         ax.annotate(
        #             vid, (updated_freq[i], panel_freq[i]), ha="center", fontsize=6
        #         )

    v_types = ["REF --> ALT", "Strand Flipped", "Ambigious Variants", "ALT Missing"]
    counts = [
        v_summary.updated,
        v_summary.flipped,
        v_summary.ambigious,
        v_summary.unknown_alt,
    ]
    bar_width = 0.75
    idx = np.arange(len(counts))
    barlist = ax3.bar(idx, counts, width=bar_width, align="center")
    ax3.set_xticks(idx)
    ax3.set_xticklabels(v_types, rotation=45, minor=False, fontsize=7)
    [bar.set_color("r") for bar in barlist[2:]]
    for i, count in enumerate(counts):
        col = "r" if i > 1 else "black"
        ax3.text(i, count, " " + str(count), ha="center", color=col)
    ax3.set_ylabel("Counts", fontsize=8)
    ax3.set_title("Variant Modification Type and Excluded Counts", fontsize=9)

    def_color = plt.rcParams["axes.prop_cycle"].by_key()["color"][0]
    leg_elements = [
        Patch(facecolor=def_color, label="Updated"),
        Patch(facecolor="red", label="Removed"),
    ]
    ax3.legend(handles=leg_elements, loc="upper right")

    (minor_freq, missing_freq, freq_diff) = zip(
        *[
            (
                v_summary.VARIANTS[v_id]["MAF"],
                v_summary.VARIANTS[v_id]["MISS"],
                v_summary.VARIANTS[v_id]["PFD"],
            )
            for v_id in v_ids
        ]
    )

    gt_coefs = (
        np.corrcoef(minor_freq, freq_diff)[1, 0],
        np.corrcoef(missing_freq, freq_diff)[1, 0],
    )

    gt_freqs = (minor_freq, missing_freq)
    gt_titles = (
        "Minor Allele Frequency Compared to Difference in Frequency Between Panel and VCF",
        "Genotype Missingness Compared to Difference in Frequency Between Panel and VCF",
    )
    y_labs = ("Minor Allele Frequency", "Missing Genotype Frequency")

    for i, ax in enumerate([ax4, ax5]):
        ax.set_title(gt_titles[i], fontsize=9)
        ax.scatter(freq_diff, gt_freqs[i], s=10, alpha=0.7)
        ax.set_ylabel(y_labs[i], fontsize=8)
        ax.set_xlabel("Panel vs VCF Frequency Difference", fontsize=8)
        ax.plot(ax.get_xlim(), ax.get_ylim(), ls="--", c=".3")
        ax.annotate(
            "corr = %.2f" % gt_coefs[i],
            (max(freq_diff) - X_OFFSET, max(gt_freqs[i]) - Y_OFFSET),
            ha="center",
        )
    plt.savefig(outfile)


def generate_panel_data(
    panel_file, chr=None, annotation="ensembl", panel_type="kirimp", separator=None
):
    f = open(panel_file, "r")
    header_line = next(f).strip()
    sep = get_separator(header_line, separator)
    header_line = [cell.replace('"', "") for cell in header_line.split(sep)]
    if panel_type == "kirimp":
        chromosome = CHROMOSOME_19_ANNOTATION[annotation]
        if header_line != KIRIMP_HEADER:
            raise TypeError(
                "If input panel type is kirimp, the panel needs to contain a comma-separated header:\n%s"
                % ",".join(KIRIMP_HEADER)
            )
    else:
        if header_line != CUSTOM_HEADER:
            raise TypeError(
                "If input panel type is custom, the panel needs to contain a comma-separated header:\n%s"
                % ",".join(CUSTOM_HEADER)
            )
    snp_dict = {
        chromosome + "_" + cells[1]
        if panel_type == "kirimp"
        else cells[0]
        + "_"
        + cells[1]: {"A0": cells[2], "A1": cells[3], "freq": float(cells[4].strip())}
        if panel_type == "kirimp"
        else {"A0": cells[2], "A1": cells[3], "freq": float(cells[4])}
        for cells in [line.strip().replace('"', "").split(sep) for line in f]
    }
    f.close()
    return snp_dict


def get_separator(line, passed_separator=None):
    tabs = line.count(r"\t")
    commas = line.count(r",")
    if passed_separator:
        sep = passed_separator
    elif tabs == 4 and commas != 4:
        sep = r"\t"
    elif tabs != 4 and commas == 4:
        sep = ","
    else:
        raise TypeError(
            "Cannot determine separator from file please specify separator directly as an argument [--reference-panel-col-separator]"
        )
    return sep


def parse_arguments(arguments=None):
    parser = argparse.ArgumentParser(
        description="This script encodes SNPs in a VCF to a reference panel based on allele frequencies",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser._action_groups.pop()
    required = parser.add_argument_group("required arguments")
    optional = parser.add_argument_group("optional arguments")
    required.add_argument(
        "-v",
        "--vcf-file",
        help="VCF/BCF file to re-encode (can be compressed with bgzip)",
        required=True,
        type=str,
    )
    required.add_argument(
        "-r",
        "--reference-panel",
        help="Reference panel file  either in format of KIR*IMP reference panel  or a custom data format [chrom pos ref alt freq]",
        required=True,
        type=str,
    )
    optional.add_argument(
        "-rt",
        "--reference-panel-type",
        help="Reference panel file type",
        choices=["kirimp", "custom"],
        default="kirimp",
        type=str,
    )
    optional.add_argument(
        "--separator", help="Custom reference panel column separator", type=str
    )
    optional.add_argument(
        "-o", "--output", help="Output vcf file", required=False, type=str
    )
    optional.add_argument(
        "-chr",
        "--chromosomes",
        help="Chromosome over which to encode SNPs ",
        required=False,
        nargs="?",
        type=str,
    )
    optional.add_argument(
        "--chromosome-annotation",
        help="Chromosome annotation type in the VCF",
        choices=["ucsc", "ensembl", "genbank"],
        default="ensembl",
        type=str,
    )
    optional.add_argument(
        "-a",
        "--ambigious",
        help="Determines whether ambigious alternate alleles should be dropped",
        action="store_false",
    )
    optional.add_argument(
        "-c",
        "--fix-complement-ref-alt",
        help="Should ref/alt that are complements be fixed with respect to frequency",
        action="store_false",
    )
    optional.add_argument(
        "-min",
        "--min-ambigious-threshold",
        help="Alternate alleles above this frequency and below the max ambigious frequency will be flagged as ambigious",
        default=0.495,
        type=float,
    )
    optional.add_argument(
        "-max",
        "--max-ambigious-threshold",
        help="Alternate alleles above this frequency and below the max ambigious frequency will be flagged as ambigious",
        default=0.505,
        type=float,
    )
    args = vars(parser.parse_args())
    if (
        args["reference_panel_type"] == "custom"
        and args["reference_panel_format"] is None
    ):
        parser.error(
            "custom --reference-panel-type requires --reference-panel-format to be set"
        )
    return args


if __name__ == "__main__":
    main()
