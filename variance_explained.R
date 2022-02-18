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
      data.frame(metabolite = x, r.squared = r.squared, n = nrow(tax), mgs = paste(res$name, collapse = ";"), 
                 importance = paste(res$coefficient, collapse = ";"), error = NA)
      
    }, error = function(e) {
      
      # return empty result if error
      
      data.frame(metabolite = NA, r.squared = NA, n = NA, mgs = NA, importance = NA, error = paste("Error:", e$message))
      
    })

  }, BPPARAM = MulticoreParam(cores))
  
  do.call(rbind, res)

}

# import data

load("data.rda")

# run lasso

variance.explained <- var.fun()

# order

variance.explained <- variance.explained[order(-abs(variance.explained$r.squared)), ]

# export data

export(variance.explained[, c("metabolite", "r.squared", "n", "mgs", "importance")],
       "variance_explained.tsv")
