# set options and load libraries

cores <- 16

library(rio)
library(BiocParallel)
library(fgsea)

module.fun <- function(data) {
  
  pathways <- module
  
  res <- bplapply(unique(data$metabolite), function(metabolite) {
    
    data <- data[which(data$metabolite == metabolite), ]
    stats <- rank(-data$p.value, na = "keep")
    names(stats) <- data$mgs
    res <- as.data.frame(fgsea(pathways, stats, scoreType = "pos", eps = 0))
    data.frame(metabolite = metabolite, module = res$pathway, estimate = res$NES, p.value = res$pval, 
               size = res$size, leading = sapply(res$leadingEdge, function(x) paste(x, collapse = ";")))
    
  }, BPPARAM = MulticoreParam(cores))
  
  do.call(rbind, res)

}

genus.fun <- function(data) {
  
  pathways <- genus
  
  res <- bplapply(unique(data$metabolite), function(metabolite) {
    
    data <- data[which(data$metabolite == metabolite), ]
    stats <- rank(-data$p.value, na = "keep")
    names(stats) <- data$mgs
    res <- as.data.frame(fgsea(pathways, stats, scoreType = "pos", eps = 0))
    data.frame(metabolite = metabolite, genus = res$pathway, estimate = res$NES, p.value = res$pval, 
               size = res$size, leading = sapply(res$leadingEdge, function(x) paste(x, collapse = ";")))
    
  }, BPPARAM = MulticoreParam(cores))
  
  do.call(rbind, res)

}

subclass.fun <- function(data) {
  
  pathways <- subclass
  
  res <- bplapply(unique(data$mgs), function(mgs) {
    
    data <- data[which(data$mgs == mgs), ]
    stats <- rank(-data$p.value, na = "keep")
    names(stats) <- data$metabolite
    res <- as.data.frame(fgsea(pathways, stats, scoreType = "pos", eps = 0))
    data.frame(mgs = mgs, subclass = res$pathway, estimate = res$NES, p.value = res$pval, 
               size = res$size, leading = sapply(res$leadingEdge, function(x) paste(x, collapse = ";")))
    
  }, BPPARAM = MulticoreParam(cores))
  
  do.call(rbind, res)

}

# import data

load("data/data.rda")
correlations <- import("results/correlations.tsv")
correlations <- correlations[which(!is.na(correlations$p.value)), ]

# test data

# correlations <- data.frame(mgs = rep(paste0("mgs", 1:30), each = 30), 
#                            metabolite = rep(paste0("metabolite", 1:30), 30), 
#                            estimate = runif(900, -1, 1), 
#                            p.value = runif(900)
#                            )
# 
# module <- list(module1 = sample(unique(correlations$mgs), 5), 
#                module2 = sample(unique(correlations$mgs), 5), 
#                module3 = sample(unique(correlations$mgs), 5)
#                )
# 
# genus <- list(genus1 = sample(unique(correlations$mgs), 5), 
#               genus2 = sample(unique(correlations$mgs), 5), 
#               genus3 = sample(unique(correlations$mgs), 5)
#               )
# 
# subclass <- list(subclass1 = sample(unique(correlations$metabolite), 5), 
#                  subclass2 = sample(unique(correlations$metabolite), 5), 
#                  subclass3 = sample(unique(correlations$metabolite), 5)
#                  )

# run gsea for positive correlations

module.pos <- module.fun(correlations[which(correlations$estimate >= 0), ])
genus.pos <- genus.fun(correlations[which(correlations$estimate >= 0), ])
subclass.pos <- subclass.fun(correlations[which(correlations$estimate >= 0), ])

# run gsea for negative correlations

module.neg <- module.fun(correlations[which(correlations$estimate < 0), ])
genus.neg <- genus.fun(correlations[which(correlations$estimate < 0), ])
subclass.neg <- subclass.fun(correlations[which(correlations$estimate < 0), ])

# add direction

module.pos$direction <- "positive"
genus.pos$direction <- "positive"
subclass.pos$direction <- "positive"

module.neg$direction <- "negative"
genus.neg$direction <- "negative"
subclass.neg$direction <- "negative"

# combine

module.res <- rbind(module.pos, module.neg)
genus.res <- rbind(genus.pos, genus.neg)
subclass.res <- rbind(subclass.pos, subclass.neg)

# adjust for multiple testing

module.res$q.value <- p.adjust(module.res$p.value, method = "fdr", n = sum(!is.na(module.res$p.value)))
genus.res$q.value <- p.adjust(genus.res$p.value, method = "fdr", n = sum(!is.na(genus.res$p.value)))
subclass.res$q.value <- p.adjust(subclass.res$p.value, method = "fdr", n = sum(!is.na(subclass.res$p.value)))

# order

module.res <- module.res[order(module.res$p.value, -abs(module.res$estimate)), ]
genus.res <- genus.res[order(genus.res$p.value, -abs(genus.res$estimate)), ]
subclass.res <- subclass.res[order(subclass.res$p.value, -abs(subclass.res$estimate)), ]

# export data

export(module.res[, c("metabolite", "module", "estimate", "p.value", "q.value", "leading", "size")], 
       "results/gsea_module.tsv")
export(genus.res[, c("metabolite", "genus", "estimate", "p.value", "q.value", "leading", "size")], 
       "results/gsea_genus.tsv")
export(subclass.res[, c("mgs", "subclass", "estimate", "p.value", "q.value", "leading", "size")], 
       "results/gsea_subclass.tsv")
