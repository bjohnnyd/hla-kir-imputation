##Steps to convert Traherne PLINK BED from b36  to VCF file of build NCBIv37/hg19

##TOOLS:

1. plink (version 1.90b4)
2. ucsc-liftover (version 377)

##DATA

1. Get `liftOver` chain file at:
  - `http://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz`

##STEPS:

1. Created `UCSC BED` file from `PLINK BED` file to be used with `liftOver`: 
  - `plink --bfile traherne --allow-no-sex --real-ref-alleles --recode --out output/traherne`
  - `awk 'BEGIN{OFS="\t"} {print "chr"$1, $4, $4 + 1, $2, $3}' output/traherne.map > output/traherne.bedfile `

2. Perform `liftOver`: 
  - `liftOver /output/traherne.bedfile hg18ToHg19.over.chain.gz output/liftover/traherne_hg19.bedfile output/liftover/traherne_hg18.unlifted`

3. Create PLINK update file:
  - `awk '{print $4, $2}' OFS='\t' output/liftover/traherne_hg19.bedfile > /output/liftover/traherne_hg19.changefile`

4. Update PLINK files to b37:
  - `plink --file output/traherne --allow-no-sex --real-ref-alleles --make-bed --update-map /output/liftover/traherne_hg19.changefile --exclude output/liftover/traherne_hg18.unlifted --out output/traherne_b37`

5. Convet PLINK file to VCF:
  - `plink --allow-extra-chr --allow-no-sex --real-ref-alleles --keep-allele-order --snps-only just-acgt --bfile output/traherne_b37 --recode vcf --out output/vcf/traherne_b37`

##FOR KIR*IMP 

The VCF file will still have a differences between the `KIR*IMP` reference and alt alleles due to differences in strand and also REF/ALT being coded opposite.  The REF/ALT in the VCF file therefore need to be first synchronized by equating REF/ALT with the `KIR*IMP` panel in cases where the REF/ALT are opposite. In cases, where the KIR*IMP panel is on opposite strand the VCF file should be flipped.  You can compare VCF to [dbSNP](https://www.ncbi.nlm.nih.gov/snp/) for genome build b37 by position to convert to forward/reverse coding.   
