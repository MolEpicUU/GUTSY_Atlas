# load libraries and set options

library(rio)
library(BiocParallel)
library(ppcor)

cores <- 16

# correlation function

cor.fun <- function(y, x, z) {
  
  fit <- apply(y, 2, function(y) {
    
    tryCatch({
      
      # create model matrix
      
      data <- data.frame(y, x, z)
      data <- data[complete.cases(data), ]
      data <- data[, apply(data, 2, function(x) length(unique(x)) > 1)]
      data <- as.data.frame(model.matrix(~ ., data))[, -1]
      
      # run model
      
      data.frame(pcor.test(data$y, data$x, data[, which(!colnames(data) %in% c("y", "x"))], method = "spearman"), 
                 message = NA)
      
    }, warning = function(w) {
      
      # if warning, run model again to capture warning

      data <- data.frame(y, x, z)
      data <- data[complete.cases(data), ]
      data <- data[, apply(data, 2, function(x) length(unique(x)) > 1)]
      data <- as.data.frame(model.matrix(~ ., data))[, -1]
      data.frame(pcor.test(data$y, data$x, data[, which(!colnames(data) %in% c("y", "x"))], method = "spearman"), 
                 message = paste("Warning:", w$message))
      
    }, error = function(e) {
      
      # return empty result if error
      
      data.frame(estimate = NA, p.value = NA, statistic = NA, n = NA, gp = NA, Method = NA, 
                 message = paste("Error:", e$message))
      
    })
    
  })
  
  fit <- do.call(rbind, fit)
  fit[, c("estimate", "p.value", "statistic", "n", "message")]
  
}

# main function

main.fun <- function() {
  
  # calculate correlations
  
  res <- bplapply(colnames(tax), function(x) {
    
    coef <- cor.fun(met, tax[, x], pheno[, c("age", "sex", "ethnicity", "batch")])
    data.frame(mgs = x, metabolite = colnames(met), coef)
    
  }, BPPARAM = MulticoreParam(cores))
  
  res <- do.call(rbind, res)
  res[, c("mgs", "metabolite", "estimate", "statistic", "p.value", "n", "message")]
  
}

# import data

load("data.rda")

# add alpha diversity

tax$shannon <- pheno$shannon

# run models

res <- main.fun()
diversity <- res[which(res$mgs == "shannon"), ]
res <- res[which(res$mgs != "shannon"), ]

# adjust for multiple testing

res$q.value <- p.adjust(res$p.value, method = "fdr", n = sum(!is.na(res$p.value)))
diversity$q.value <- p.adjust(diversity$p.value, method = "fdr", n = sum(!is.na(diversity$p.value)))

# order

res <- res[order(res$p.value, -abs(res$estimate)), ]
diversity <- diversity[order(diversity$p.value, -abs(diversity$estimate)), ]

# export data

export(res[, c("mgs", "metabolite", "estimate", "statistic", "p.value", "q.value", "n", "message")], 
       "correlations.tsv")
export(diversity[, c("metabolite", "estimate", "statistic", "p.value", "q.value", "n", "message")], 
       "diversity.tsv")
