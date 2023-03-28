library(MASS)
library(corrplot)
library(Hmisc)
library(factoextra)

locCalc <- function(dataset){
	loc <- c()
	for(i in 1:length(dataset[,1])){
		gene <- subset(genes, name==dataset[i, ]$gene)
		if(dataset[i, ]$effect %in% c("upstream", "UTR-5")){
			if(gene$Strand == "+") loc[i] <- paste0("+", gene$chromStart - dataset[i, ]$pos)
			else loc[i] <- paste0("+", dataset[i, ]$pos - gene$chromEnd)
		}
	if(dataset[i, ]$effect %in% c("downstream", "UTR-3")){
    if(gene$Strand == "+") loc[i] <- paste0("-", dataset[i,]$pos - gene$chromEnd)
		else loc[i] <- paste0("-", gene$chromStart - dataset[i, ]$pos)
	}
	if(dataset[i, ]$effect %nin% c("upstream", "UTR-5", "UTR-3", "downstream")){
    if(gene$Strand == "+") loc[i] <- dataset[i,]$pos - gene$chromStart
    else loc[i] <- gene$chromEnd - dataset[i,]$pos
}
}
	loc
}


##########################
#Load strain information
##########################
strains <- read.csv("data_in/general/180727strainInfo.csv")
strains <- subset(strains, ST=="93" & arm != "pre-COAT")

centelo <- read.table("data_in/general/centromeres_telomeres_cna3.gff", sep="\t", as.is=TRUE)
centelo <- rbind(centelo, c("chr5", "prediction", "MAT", 185700 , 275730, ".", "+", ".", "locus_tag Chr5 MAT"))

chrLens <- c(2291499, 1621675, 1575141, 1084805, 1814975, 1422463, 1399503, 1398693, 1186808, 1059964, 1561994, 774062, 756744, 926563)

genes <- read.csv("data_in/general/180709genesList.csv")

###########################
#Read in phenotype files
###########################
dCSF <- read.csv("data_in/phenotypes/Clean_Compiled_CSF_log2duplicates_d0.csv", as.is=TRUE) #"duplicates", t=0

dWBC <- read.csv("data_in/phenotypes/COAT_WBC.csv")
dWBC.sub <- subset(dWBC, line %in% dCSF$line)

invitro.sub <- subset(strains, line %in% dCSF$line )

#exclude LFA titer, protein, VEGF, HIV-viral and CSF WBC
df <- data.frame(dCSF[, c(1:21)], dWBC.sub[, c(5, 7, 9)], invitro.sub[, c(5:11)])

df.sub <- na.omit(df) #29 lines remain

###########################
#Read in variant files
###########################
var.dCSF <- read.csv("data_out/variants/var_dCSF.csv")
var.dCSF$CP <- paste(var.dCSF$CHROM, var.dCSF$POS, sep=".")
var.dWBC <- read.csv("data_out/variants/var_dWBC.csv")
var.dWBC$CP <- paste(var.dWBC$CHROM, var.dWBC$POS, sep=".")

###########################
#PCA analysis
###########################
# singular value decomposition of the centered and scaled to have unit variance data
#remove CSF-WBC and the MIC tests (missing values mean only 24 lines could be used


prin_comp <- prcomp(df.sub[, 5:31], scale. = T)
#prin_comp <- prcomp(~V1+V2, data=d, center = TRUE, scale = TRUE, na.action = na.omit)
biplot(prin_comp, scale = 1, pc.biplot=TRUE)

#x = the value of the rotated data (the centred (and scaled if requested) data multiplied by the rotation matrix) is returned. Hence, cov(x) is the diagonal matrix diag(sdev^2)
PCAvar <-data.frame(df.sub[,1:4], prin_comp$x)
write.csv(PCAvar, "data_out/PCA/180709PCAloadings.csv", row.names=FALSE)

PCAvar$ST93 <- factor(PCAvar$ST93, levels=c("D", "C", "A", "B"))

pdf("manuscript/figures/FigureS2B_PC12-clade.pdf", width=3.5, height=4)
plot(as.numeric(PCAvar$ST93), PCAvar$PC1, xaxt= "n", yaxt="n", xlab="ST93 clade", ylim=c(-8, 8), col="red", ylab = "PC value")
points(as.numeric(PCAvar$ST93)+0.2, PCAvar$PC2, col="blue")
axis(1, 1:4, c("UgCl362", "UgCl495", "Clade A", "Clade B"), cex.axis=0.8)
axis(2, las=2, at=c(-8, -4, 0, 4, 8), cex.axis=0.9)
legend("topleft", col=c("red", "blue"), legend = c("PC1", "PC2"), pch=21)
dev.off()


#compute standard deviation of each principal component
std_dev <- prin_comp$sdev

#compute variance (eigenvalues)
pr_var <- std_dev^2

prop_varex <- pr_var/sum(pr_var)

plot(pr_var, xlab = "Principal Component",
            ylab = "Eigenvalues",
            type = "b", pch=19, yaxt="n", xaxt="n", cex=1.2)
 axis(2, las=2, cex.axis=0.8)
axis(1, at=1:21, cex.axis=0.8)

#http://www.quantpsy.org/pubs/preacher_maccallum_2003.pdf
df.rand <- apply(df.sub[,5:31], 2, sample)
prin_comp.rand <- prcomp(df.rand, scale. = T)
PCAvar.rand <-data.frame(df.sub[,1:5], prin_comp.rand$x)
std_dev.rand <- prin_comp.rand$sdev
pr_var.rand <- std_dev.rand^2
prop_varex.rand <- pr_var.rand/sum(pr_var.rand)

points(pr_var.rand, pch=21, type="b", lty=2, cex=1.2)
legend("topright", pch=c(19, 21), legend=c("actual data", "randomized data"), lty=1)

pdf("manuscript/figures/FigureS2A_propVarex-data-vs-20rand.pdf", width=6, height=4)
plot(prop_varex*100, xlab = "Principal Component",
            ylab = "Percent of Variance Explained",
            type = "b", pch=19, yaxt="n", xaxt="n", cex=1.2, col="red")
 axis(2, las=2, cex.axis=0.7)
axis(1, at=seq(1, 27, by=2), cex.axis=0.8)
axis(1, at=seq(2, 26, by=2), cex.axis=0.8)

#http://www.quantpsy.org/pubs/preacher_maccallum_2003.pdf
i <- 0
while(i <11){
  i <- i+1
  df.rand <- apply(df.sub[,5:31], 2, sample)
  prin_comp.rand <- prcomp(df.rand, scale. = T)
  std_dev.rand <- prin_comp.rand$sdev
  pr_var.rand <- std_dev.rand^2
  prop_varex.rand <- pr_var.rand/sum(pr_var.rand)
  points(prop_varex.rand*100, pch=21, type="l", lty=2, cex=0.8)
}
points(prop_varex*100, col="red", type="b", cex=1.2, pch=19)
legend("topright", pch=c(19, 0), legend=c("actual data", "randomized data"), lty=1)
dev.off()

######################################################################
#use PCs in logistic regression calculation
######################################################################
#var.dCSF.sub <- var.dCSF[, c(2:12, 14:20, 22:34, 36:48)]
#exclude 8 strains (UgCl212, UgCl332, UgCl357, UgCl422, UgCl447, UgCl461, UgCl541, UgCl549
#same as df.sub$line
var.dCSF.sub <- var.dCSF[, c(2:12, 15:20, 22, 24:34, 36:48)]

j <- 0
pvalPCA <- list()
for(i in c(5:6)){ #PC1 and PC2
	j <- j+1
	pvals.PCA <- c()
	for(k in 1:length(var.dCSF.sub[,1])){
	      resp <- t(var.dCSF.sub[k,1:29])
	      pred <- t(PCAvar[,i])
				dat <- data.frame(pred=t(pred), resp)
	      names(dat) <- c("pred", "resp")
	      dat <- subset(dat, !is.na(pred))
	      dat <- subset(dat, !is.na(resp))
	      dat$resp <- as.character(dat$resp)
	      dat$resp[dat$resp!= var.dCSF.sub[k,32]] <- 1		#ALT
	      dat$resp[dat$resp== var.dCSF.sub[k,32]] <- 0		#REF
	      dat$pred <- as.numeric(dat$pred)
	      dat$resp <- as.numeric(dat$resp)
	      t2 <- glm(resp~pred, data=dat, family=binomial)
				pvals.PCA[k] <- summary(t2)$coefficients[2,4]
	}
	pvalPCA[[j]] <- pvals.PCA
}

names(pvalPCA) <- c("PCA1", "PCA2")


#geneP.PCA <- data.frame(var.dCSF.sub[,32:44], "PCA1" = unlist(pvalPCA[[1]]), "PCA2" = unlist(pvalPCA[[2]]), "PCA3" = unlist(pvalPCA[[3]]))
geneP.PCA <- data.frame(var.dCSF.sub[, 30:42], "PCA1" = unlist(pvalPCA[[1]]), "PCA2" = unlist(pvalPCA[[2]]))
names(geneP.PCA)[1:12] <- c("chr", "pos", "REF", "ALT", "QUAL", "HOM.REF", "HOM.VAR", "gene", "effect", "impact", "class", "aa_change")
geneP.PCA.sig <- subset(geneP.PCA, PCA1 < 0.05 | PCA2 < 0.05 )

geneP.PCA.sig$gene <- as.character(geneP.PCA.sig$gene)

alias <- c()
description <- c()
for(i in 1:length(geneP.PCA.sig$gene)){
  if(is.na(as.character(geneP.PCA.sig$gene[i]))) alias[i] <- "null"
  else alias[i] <- as.character(subset(genes, name==as.character(geneP.PCA.sig$gene[i]))$alias)
  description[i] <- as.character(subset(genes, name==as.character(geneP.PCA.sig$gene[i]))$description)
}

geneP.PCA.sig$alias <- alias
geneP.PCA.sig$description <- description
geneP.PCA.sig <- geneP.PCA.sig[order(geneP.PCA.sig$gene),]

geneP.PCA.sig$effect <- as.character(geneP.PCA.sig$effect)
geneP.PCA.sig$effect[geneP.PCA.sig$effect=="NON_SYNONYMOUS_CODING"] <- "ns"
geneP.PCA.sig$effect[geneP.PCA.sig$effect=="START_GAINED"] <- "start+"
geneP.PCA.sig$effect[geneP.PCA.sig$effect=="UTR_5_PRIME"] <- "UTR-5"
geneP.PCA.sig$effect[geneP.PCA.sig$effect=="UTR_3_PRIME"] <- "UTR-3"
geneP.PCA.sig$effect[geneP.PCA.sig$effect=="UPSTREAM"] <- "upstream"
geneP.PCA.sig$effect[geneP.PCA.sig$effect=="FRAME_SHIFT"] <- "frameshift"
geneP.PCA.sig$effect[geneP.PCA.sig$effect=="DOWNSTREAM"] <- "downstream"

geneP.PCA.sig$gene <- as.character(geneP.PCA.sig)

#Table3
write.csv(geneP.PCA.sig, "manuscript/tables/Table3_180708geneP-PCA_sig.csv", row.names=FALSE)
