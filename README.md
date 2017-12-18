
<!-- README.md is generated from README.Rmd. Please edit that .Rmd file -->
pKSEA
=====

The goal of pKSEA is to infer kinase activity from phosphoproteomics data using in-silico kinase-substrate predictions. pKSEA uses summary statistics calculated from phosphoproteomic data at the peptide level to infer changes in kinase activity across experimental conditions. pKSEA then uses kinase-substrate prediction scores to weight observed changes in phosphopeptide abundance to calculate a phosphopeptide-level contribution score, then sums up these contribution scores by kinase to obtain a phosphoproteome-level kinase activity change score (KAC score). pKSEA then assesses the significance of changes in predicted substrate abundances for each kinase using permutation testing. This results in a permutation score (pKSEA significance score) reflecting the likelihood of a similarly high or low KAC from random chance, which can then be interpreted in an analogous manner to an empirically calculated p-value. pKSEA contains default databases of kinase-substrate predictions from NetworKIN (NetworKINPred\_db) and of known kinase-substrate links from PhosphoSitePlus (KSEAdb).

Please see package details and individual function information for input data formatting and additional examples.

Installation
------------

You can install pKSEA from github with:

``` r
# install.packages("devtools")
devtools::install_github("pll21/pKSEA")
```

References
----------

Horn et al., KinomeXplorer: an integrated platform for kinome biology studies. Nature Methods 2014 Jun;11(6):603–4.

Hornbeck PV, Zhang B, Murray B, Kornhauser JM, Latham V, Skrzypek E PhosphoSitePlus, 2014: mutations, PTMs and recalibrations. Nucleic Acids Res. 2015 43:D512-20.

<!-- ## Example -->
<!-- This is a basic example which shows you how to solve a common problem: -->
<!-- ```{r example} -->
<!-- #Load relevant databases into current environment -->
<!-- data("KSEAdb") -->
<!-- data("NetworKINPred_db") -->
<!-- ``` -->
