# set options and load libraries

setwd("//argos.rudbeck.uu.se/MyGroups$/Gold/MolEpi/Personal folders and scripts/Koen/manuscripts/manuscript1/scripts/demo")
set.seed(1)

library(scales)

# create test data

tax <- as.data.frame(matrix(runif(1000 * 45), ncol = 45))
colnames(tax) <- paste0("mgs", 1:ncol(tax))
rownames(tax) <- paste0("sample", 1:nrow(tax))

met <- tax + rnorm(1000 * 45) / 2
colnames(met) <- paste0("metabolite", 1:ncol(met))
rownames(met) <- paste0("sample", 1:nrow(met))
met <- apply(met, 2, function(x) {
  
  sample <- sample(1:2, 1)
  scale <- list(rescale(x, c(0, 5)), rescale(x, c(5, 0)))
  scale[[sample]]
  
})

tax2 <- lapply(1:15, function(i) {
  
  rescale(met[, 1] + rnorm(1000), c(0, 1))

})
tax2 <- do.call(cbind, tax2)
tax[, 1:15] <- tax2

pheno <- data.frame(age = rnorm(1000, 55, 4),
                    sex = sample(c("female", "male"), 1000, replace = T),
                    ethnicity = sample(c("country1", "country2", "country3"), 1000, replace = T),
                    site = sample(c("site1", "site2"), 1000, replace = T),
                    batch = sample(c("batch1", "batch2", "batch3"), 1000, replace = T),
                    shannon = abs(rnorm(1000, 3, 2))
                    )
rownames(pheno) <- paste0("sample", 1:nrow(tax))

met1 <- met[, 1:40]
met2 <- met[, 41:45]
met2[which(pheno$site == "site1"), ] <- NA
met <- cbind(met1, met2)

modules <- list(module1 = colnames(tax)[1:15],
               module2 = colnames(tax)[16:30],
               module3 = colnames(tax)[31:45]
               )

genera <- list(genus1 = colnames(tax)[1:15],
              genus2 = colnames(tax)[16:30],
              genus3 = colnames(tax)[31:45]
              )

subclasses <- list(subclass1 = colnames(met)[1:15],
                 subclass2 = colnames(met)[16:30],
                 subclass3 = colnames(met)[31:45]
                 )

# export data

save(tax, met, pheno, modules, genera, subclasses, file = "data.rda")
