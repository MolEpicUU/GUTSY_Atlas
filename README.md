Analysis code for "An online atlas of human plasma metabolite signatures of gut microbiome composition" by Dekkers et al.

Scripts:
* correlations.R is the R-script used to calculate partial Spearman correlations between gut microbial species and plasma metabolites.
* gsea.R is the R-script used to generate the enrichment test results for gut metabolic modules, taxonomic genera and metabolite classes.
* variance_explained.R is the R-script used to calculate how much of the variance of gut microbial species is explained by a combination of metabolites.

Software versions tested:
R 4.1.1
rio 0.5.27
BiocParallel 1.26.2
ppcor 1.1
glmnet 4.1-3
fgsea 1.19.2

NOTE: the scripts contain a section called "test data" that can be uncommented to test the scripts without access to the data.
