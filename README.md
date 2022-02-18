# Analysis code for "An online atlas of human plasma metabolite signatures of gut microbiome composition" by Dekkers et al.

# Scripts

* correlations.R is the R-script used to calculate partial Spearman correlations between gut microbial species and plasma metabolites
* gsea.R is the R-script used to generate the enrichment test results for gut metabolic modules, taxonomic genera and metabolite classes
* variance_explained.R is the R-script used to calculate how much of the variance of gut microbial species is explained by a combination of metabolites
* data.rda is simulated test data, the scripts work both on the real data as well as the test data

# Installation guide

* Install R: https://www.r-project.org/
* Install Bioconductor: https://bioconductor.org/install/
* Install packages in R:
  - library(BiocManager)
  - install(c("rio", "BiocParallel", "ppcor", "glmnet", "fgsea"))
* Download correlations.R, gsea.R, variance_explained.R

* Install time: <10 minutes

# Demo

* Download data.rda
* In R:
  - Set working directory to directory with data.rda, e.g. setwd("C:/Users/User/Documents/Demo/")
  - Run script correlations.R
  - Run script gsea.R
  - Run script variance_explained.R
* Expected output: 
  - Files: correlations.tsv, diversity.tsv, gsea_module.tsv, gsea_genus.tsv, gsea_subclass.tsv and variance_explained.tsv
  - In correlations.tsv mgs 1 to 30 should be associated with metabolite 1 to 30, respectively
  - In gsea_module.tsv and gsea_genus.tsv the strongest enriched feature for metabolite 1 to 30 should often contain mgs 1 to 30 as leading edge, respectively
  - In gsea_subclass.tsv the strongest enriched feature for mgs 1 to 30 should often contain metabolite 1 to 30 as leading edge, respectively
  - In variance_explained.tsv metabolite 1 to 30 should have mgs 1 to 30 as most important contributor to the variance explained, respectively
* Run time: <5 minutes

Software versions tested: R 4.1.1, rio 0.5.27, BiocParallel 1.26.2, ppcor 1.1, glmnet 4.1-3 and fgsea 1.19.2


