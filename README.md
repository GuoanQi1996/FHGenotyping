# Functional Haplotype Genotyping (FHG)
Functional Haplotype Genotyping (FHG) is design to genotype the functional haplotype (FH), represents the variations in the coding region of gene / protein, from genomic variations. Technically, assuming a sample size of _N_, a number of markers represented by SNPs and various short variants of _M_, and the number of annotated genes in the target genome of _G_, FHG converts the genotype matrix of _N_ * _M_ to the genotype matrix of _N_ * _G_. The new genotype encoding takes the gene or protein itself as the object of analysis, thus making it easier to interpret a series of subsequent genetic analyses. In particular, genome-wide association analysis of FH encoded by FHG will directly map candidate genes associated with traits.

[![Maintained - Cotton Presicion Breeding Academy](https://img.shields.io/badge/Maintained-Cotton_Presicion_Breeding_Academy-green)](http://cotton.zju.edu.cn/)

## Installation
FHG is a R script tool. Please make sure the `R` is available in the OS. Directly download the repository and the run `FHG.R` by `Rscript`.
```bash
Rscript FHG.R
```
or 
```bash
Rscript FHG.R --help
```
The packages `data.table` is required for the script, `FHG` will automatically detect the dependencies and install unavailable packages. The users can also install the required packages manually if the network is not connected.

## Input
FHG require variant-major additive component file (TRAW) as the input format, this format can be generated from the standard Variant Call Format (VCF) by [PLINK](https://www.cog-genomics.org/plink/) or [PLINK2](https://www.cog-genomics.org/plink/2.0/), by the following command. 
```bash
# PLINK 1.9
gzip -d -k test.vcf.gz
plink --vcf test.vcf --allow-extra-chr --recode A-transpose --out test
# PLINK 2.0 (recommanded, fast and easy)
plink2 --vcf test.vcf.gz --allow-extra-chr --export Av --out test
```
