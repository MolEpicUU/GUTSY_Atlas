# load libraries and set options

library(rio)
library(BiocParallel)
library(glmnet)
library(caret)

cores <- 16
set.seed(1)

# variance explained function

var.fun <- function() {
  
  res <- bplapply(colnames(met), function(x) {
    
    # create data for complete cases
    
    data <- data.frame(met = met[, x], tax)
    complete <- complete.cases(data)
    data <- data[complete, ]
    met <- data$met
    tax <- data[, which(!colnames(data) == "met")]
    
    tryCatch({
      
      # create 10 folds
      
      folds <- createFolds(met, k = 10)
      
      nested.cv <- lapply(folds, function(fold) {
        
        # create test data
        
        met_test <- met[fold]
        tax_test <- tax[fold, ]
        
        # create cross-validation data
        
        met_cv <- met[-fold]
        tax_cv <- tax[-fold, ]
        
        # create 10 cross-validation folds
        
        cv_fold <- createFolds(met_cv, k = 10)[[1]]
        
        # create validation data
        
        met_validate <- met[cv_fold]
        tax_validate <- tax[cv_fold, ]
        
        # create training data
        
        met_train <- met[-cv_fold]
        tax_train <- tax[-cv_fold, ]
        
        # fit ridge regression model
        
        fit <- glmnet(x = as.matrix(tax_train), y = met_train, alpha = 0)
        
        # determine model with minimum lambda
        
        pred <- predict(fit, newx = as.matrix(tax_validate), s = fit$lambda)
        mse <- colMeans((pred - met_validate)^2) / var(met_validate)
        lambda <- fit$lambda[which.min(mse)]
        mse.train <- min(mse)
        
        # test model with minimum lambda on test fold
        
        pred <- as.numeric(predict(fit, newx = as.matrix(tax_test), s = lambda))
        mse.test <- mean((pred - met_test)^2) / var(met_test)
        r.squared <- cor(pred, met_test)^2
        data.frame(r.squared = r.squared, lambda = lambda, mse.train = mse.train, mse.test = mse.test)
        
      })
      nested.cv <- do.call(rbind, nested.cv)
      
      data.frame(metabolite = x, r.squared = mean(nested.cv$r.squared), mse.train = mean(nested.cv$mse.train), mse.test = mean(nested.cv$mse.test), n = nrow(tax), error = NA)
      
    }, error = function(e) {
      
      # return empty result if error
      
      data.frame(metabolite = NA, r.squared = NA, mse.train = NA, mse.test = NA, n = NA, error = paste("Error:", e$message))
      
    })
    
  }, BPPARAM = MulticoreParam(cores))
  
  do.call(rbind, res)
  
}

# import data

load("data.rda")

# run nested 10-fold cross-validated ridge regression

variance.explained <- var.fun()

# order

variance.explained <- variance.explained[order(-abs(variance.explained$r.squared)), ]

# export data

export(variance.explained[, c("metabolite", "r.squared", "mse.train", "mse.test", "n")],
       "variance_explained.tsv")