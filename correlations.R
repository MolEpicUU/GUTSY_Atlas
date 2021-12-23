# set options and load libraries

cores <- 16

library(rio)
library(BiocParallel)
library(ppcor)

# correlation function

cor.fun <- function(y, x, z) {
  
  # create model matrix
  
  data <- data.frame(x, z)
  data <- as.data.frame(model.matrix(~ ., data))
  
  fit <- apply(y, 2, function(ycol) {
    
    tryCatch({
      
      # run model on complete cases
      
      data$ycol <- ycol[match(rownames(data), rownames(y))]
      data <- data[complete.cases(data), ]
      data.frame(pcor.test(data$ycol, data$x, data[, which(!colnames(data) %in% c("ycol", "x"))]), message = NA)
      
    }, warning = function(w) {
      
      # run model again if warning
      
      data$ycol <- ycol[match(rownames(data), rownames(y))]
      data <- data[complete.cases(data), ]
      data.frame(pcor.test(data$ycol, data$x, data[, which(!colnames(data) %in% c("ycol", "x"))]), 
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
  
  # rank metabolites and metagenomic species
  
  met <- apply(met, 2, function(x) rank(x, "keep"))
  tax <- apply(tax, 2, function(x) rank(x, "keep"))
  
  # handle metabolites available only in 1 site
  
  met1 <- met[, apply(met, 2, function(x) any(tapply(x, pheno$site, function(y) sum(is.na(y)) == length(y))))]
  met2 <- met[, which(!colnames(met) %in% colnames(met1))]
  
  # calculate correlations for metabolites only available in 1 site
  
  res <- bplapply(colnames(tax), function(x) {
    
    if (x == "shannon") {
      
      coef <- cor.fun(met1, tax[, x], pheno[, c("age", "sex", "ethnicity", "plate")])
      
    } else {
      
      coef <- cor.fun(met1, tax[, x], pheno[, c("age", "sex", "ethnicity", "shannon", "plate")])
      
    }
    
    data.frame(mgs = x, metabolite = colnames(met1), coef)
    
  }, BPPARAM = MulticoreParam(cores))
 
  res <- do.call(rbind, res)
  
  # calculate correlations for metabolites only available in 1 site
  
  res2 <- bplapply(colnames(tax), function(x) {
    
    if (x == "shannon") {
      
      coef <- cor.fun(met2, tax[, x], pheno[, c("age", "sex", "ethnicity", "site", "plate")])
      
    } else {
      
      coef <- cor.fun(met2, tax[, x], pheno[, c("age", "sex", "ethnicity", "site", "shannon", "plate")])
      
    }
    
    data.frame(mgs = x, metabolite = colnames(met2), coef)
    
  }, BPPARAM = MulticoreParam(cores))
  
  res2 <- do.call(rbind, res2)
  res <- rbind(res, res2)
  res[, c("mgs", "metabolite", "estimate", "statistic", "p.value", "n", "message")]
  
}

# import data

load("data/data.rda")

# test data

# tax <- as.data.frame(matrix(runif(1000 * 10), ncol = 10))
# colnames(tax) <- paste0("mgs", 1:ncol(tax))
# rownames(tax) <- paste0("sample", 1:nrow(tax))
# 
# met <- as.data.frame(matrix(rnorm(1000 * 10, 5, 2), ncol = 10))
# colnames(met) <- paste0("metabolite", 1:ncol(met))
# rownames(met) <- paste0("sample", 1:nrow(met))
# 
# pheno <- data.frame(age = rnorm(1000, 55, 4), 
#                     sex = sample(c("female", "male"), 1000, replace = T), 
#                     ethnicity = sample(c("country1", "country2", "country3"), 1000, replace = T),
#                     site = sample(c("site1", "site2"), 1000, replace = T),
#                     shannon = abs(rnorm(1000, 3, 2)),
#                     plate = sample(c("plate1", "plate2", "plate3"), 1000, replace = T)
#                     )
# rownames(pheno) <- paste0("sample", 1:nrow(tax))
# 
# met1 <- met[, 1:8]
# met2 <- met[, 9:10]
# met2[which(pheno$site == "site1"), ] <- NA
# met <- cbind(met1, met2)
# rm(met1, met2)

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
       "results/correlations.tsv")
export(diversity[, c("metabolite", "estimate", "statistic", "p.value", "q.value", "n", "message")], 
       "results/diversity.tsv")
