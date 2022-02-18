# set options and load libraries

set.seed(1)
cores <- 16

library(rio)
library(BiocParallel)
library(fgsea)

# enrichment function for gmm modules

module.fun <- function(data) {
  
  res <- bplapply(unique(data$metabolite), function(metabolite) {
    
    # rank based on p-value
    
    data <- data[which(data$metabolite == metabolite), ]
    stats <- rank(-data$p.value, na = "keep")
    names(stats) <- data$mgs
    
    # run gsea
    
    res <- as.data.frame(fgsea(module, stats, scoreType = "pos", eps = 0))
    data.frame(metabolite = metabolite, module = res$pathway, estimate = res$NES, p.value = res$pval, 
               size = res$size, leading = sapply(res$leadingEdge, function(x) paste(x, collapse = ";")))
    
  }, BPPARAM = MulticoreParam(cores))
  
  do.call(rbind, res)

}

# enrichment function for genera

genus.fun <- function(data) {
  
  res <- bplapply(unique(data$metabolite), function(metabolite) {
    
    # rank based on p-value
    
    data <- data[which(data$metabolite == metabolite), ]
    stats <- rank(-data$p.value, na = "keep")
    names(stats) <- data$mgs
    
    # run gsea
    
    res <- as.data.frame(fgsea(genus, stats, scoreType = "pos", eps = 0))
    data.frame(metabolite = metabolite, genus = res$pathway, estimate = res$NES, p.value = res$pval, 
               size = res$size, leading = sapply(res$leadingEdge, function(x) paste(x, collapse = ";")))
    
  }, BPPARAM = MulticoreParam(cores))
  
  do.call(rbind, res)

}

# enrichment function for metabolite subclasses

subclass.fun <- function(data) {
  
  res <- bplapply(unique(data$mgs), function(mgs) {
    
    # rank based on p-value
    
    data <- data[which(data$mgs == mgs), ]
    stats <- rank(-data$p.value, na = "keep")
    names(stats) <- data$metabolite
    
    # run gsea
    
    res <- as.data.frame(fgsea(subclass, stats, scoreType = "pos", eps = 0))
    data.frame(mgs = mgs, subclass = res$pathway, estimate = res$NES, p.value = res$pval, 
               size = res$size, leading = sapply(res$leadingEdge, function(x) paste(x, collapse = ";")))
    
  }, BPPARAM = MulticoreParam(cores))
  
  do.call(rbind, res)

}

# import data

load("data.rda")
correlations <- import("correlations.tsv")
correlations <- correlations[which(!is.na(correlations$p.value)), ]

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
       "gsea_module.tsv")
export(genus.res[, c("metabolite", "genus", "estimate", "p.value", "q.value", "leading", "size")], 
       "gsea_genus.tsv")
export(subclass.res[, c("mgs", "subclass", "estimate", "p.value", "q.value", "leading", "size")], 
       "gsea_subclass.tsv")
