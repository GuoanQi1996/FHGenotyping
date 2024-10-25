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
The packages `data.table` is required for the script, `FHG` will automatically detect the dependencies and install an unavailable packages. The users can also install the required packages manually if the network is not connected.
