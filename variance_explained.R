# set options and load libraries

set.seed(1)
cores <- 8

library(rio)
library(BiocParallel)
library(glmnet)

# variance explained function

var.fun <- function() {
  
  res <- bplapply(colnames(met), function(x) {
    
    # create data for complete cases
    
    complete <- complete.cases(cbind(met[, x], tax))
    tax <- tax[complete, ]
    met <- met[complete, x]
    
    tryCatch({
      
      # fit model
      
      fit <- cv.glmnet(x = as.matrix(tax), y = met, alpha = 1)
      
      # calculate r.squared from cross validated errors
      
      r.squared <- 1 - fit$cvm[which(fit$lambda == fit$lambda.min)] / var(met)
      
      # calculate normalized regression coefficients
      
      coef <- coef(fit, s = "lambda.min")
      sds <- apply(as.matrix(tax), 2, sd)
      std_coef <- coef[-1, 1] * sds
      
      # clean
      
      res <- data.frame(name = colnames(tax), coefficient = std_coef)
      res <- res[which(res$coefficient != 0), ]
      res <- res[order(abs(res$coefficient), decreasing = T), ]
      data.frame(metabolite = x, r.squared = r.squared, n = n, mgs = paste(res$name, collapse = ";"), 
                 importance = paste(res$coefficient, collapse = ";"))
      
    }, error = function(e) {
      
      # return empty result if error
      
      data.frame(metabolite = NA, r.squared = NA, n = NA, mgs = NA, importance = NA)
      
    })

  }, BPPARAM = MulticoreParam(cores))
  
  do.call(rbind, res)

}

# import data

load("data/data.RData")

# test data

# tax <- as.data.frame(matrix(runif(1000 * 10), ncol = 10))
# colnames(tax) <- paste0("mgs", 1:ncol(tax))
# rownames(tax) <- paste0("sample", 1:nrow(tax))
# 
# met <- tax + rnorm(10000) # mgs1 is correlated with metabolite1, mgs2 is correlated with metabolite2, etc.
# colnames(met) <- paste0("metabolite", 1:ncol(met))
# rownames(met) <- paste0("sample", 1:nrow(met))

# run lasso

variance.explained <- var.fun()

# order

variance.explained <- variance.explained[order(-abs(variance.explained$r.squared)), ]

# export data

export(variance.explained[, c("metabolite", "r.squared", "n", "mgs", "importance")],
       "results/variance_explained.tsv")
