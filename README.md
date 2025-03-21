# Functional Genotypic Configuration Genotyping (FGCG)
Functional Genotypic Configuration Genotyping (FGCG) is design to genotype the functional genotypic configuration (FGC), represents the variations in the coding region of gene / protein, from genomic variations. Technically, assuming a sample size of _N_, a number of markers represented by SNPs and various short variants of _M_, and the number of annotated genes in the target genome of _G_, FGCG converts the genotype matrix of _N_ * _M_ to the genotype matrix of _N_ * _G_. The new genotype encoding takes the gene or protein itself as the object of analysis, thus making it easier to interpret a series of subsequent genetic analyses. In particular, genome-wide association analysis of FGC encoded by FGCG will directly map candidate genes associated with traits.

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
FHG requires variant-major additive component file (TRAW) as the input format, this format can be generated from the standard Variant Call Format (VCF) by [PLINK](https://www.cog-genomics.org/plink/) or [PLINK2](https://www.cog-genomics.org/plink/2.0/), by the following command. 
```bash
# PLINK 1.9
gzip -d -k test.vcf.gz
plink --vcf test.vcf --allow-extra-chr --recode A-transpose --out test
# PLINK 2.0 (recommanded, fast and easy)
plink2 --vcf test.vcf.gz --allow-extra-chr --export Av --out test
```
FHG also requires three annotation files for parameters `--ann`, `--mis` and `--syn`, which specify the variant annotation and defined nonsynonymous/synonymous mutation types. For example, for variant annotation information, a three-columns text file is required:
```bash
SNP	Gene	Type
A01_1945_G_A	GH_A01G0001	missense_variant
A01_2006_A_G	GH_A01G0001	synonymous_variant
A01_2010_G_A	GH_A01G0001	missense_variant
A01_2011_C_A	GH_A01G0001	missense_variant
A01_2128_T_A	GH_A01G0001	stop_gained
A01_17767_C_T	GH_A01G0002	synonymous_variant
...
```
For file defines the type of nonsynonymous mutation, a one column text file is required:
```bash
disruptive_inframe_deletion
disruptive_inframe_insertion
frameshift_variant
missense_variant
start_lost
...
```
And similarly, for file defines the type of synonymous mutation, a one column text file is required:
```bash
synonymous_variant
...
```
## Running
After installing the required runtime environment and preparing the required files, the functional genetic configuration can be generated by the following command (find the exampled dataset in `Examples` folder):
```bash
Rscript FHG.R --tfile test.traw --ann Variants_annotation.txt --mis Non_synonymous_annotation.txt --syn Synonymous_annotation.txt --out Test
```
As you may awared, this genotyping strategy can be easily extended to any genomic regions of interests, by providing proper annotation files to the software. For example, if we focus on the genomic regions located the the upsteam 2K bp region, we can annotate these variants to like `regulatory region variant` in file provided to parameter `--ann` and define the `regulatory region variant` as coding object in file provided to parameter `--mis`.
