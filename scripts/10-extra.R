library(tidyverse)

#https://www.datacamp.com/community/tutorials/survival-analysis-R

#This is the GWAS results
t <- read.csv("manuscript/tables/TableS6_geneP-allSigCP.csv")

t.ddn <- split(t, t$gene)

byGene <- data.frame(gene= names(t.ddn))
byGene$gene <- as.character(byGene$gene)
byGene$numVar <- unlist(lapply(t.ddn, function(x) nrow(x)))
byGene$numPhen <- unlist(lapply(t.ddn, function(x) sum(x$numPhen)))

#This is all sig genes (GWAS + PCA)
t2 <- read_csv("manuscript/tables/TableS9_UgCl-allSigThings.csv")
t2.gwas <- subset(t2, Gene %in% byGene$gene)
names(t2.gwas)[6] <- "KOSurv"

byGene <- cbind(byGene, t2.gwas)

tVir <- subset(byGene,  !is.na(byGene$KOSurv))
tVir$KONum <- ifelse(tVir$KOSurv == "no effect", 0, 1)

test <- glm(KONum~numVar, data=tVir, family=binomial)
tidy(test)
