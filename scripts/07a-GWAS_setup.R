library(Hmisc)
library(stringr)
library(stringi)
library(epiDisplay)

upstreamCalc <- function(dataset){
	up <- c()
	for(i in 1:length(dataset[,1])){
		gene <- subset(genes, name==dataset[i, ]$gene)
		if(dataset[i, ]$POS > gene$chromEnd) up[i] <- dataset[i, ]$POS-gene$chromEnd
		else up[i] <- gene$chromStart-dataset[i, ]$POS
	}
	up
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
#there are 35 lines with CSF data: missing 447, 461, 541 (+230, 393, 443)
dCSF$IL1b <- as.numeric(dCSF$IL1b)
# dCSF93.dup.0.sub <- subset(dCSF93.dup.0, line!="468")
#Not enough data for MIP1a, VEGF

dWBC <- read.csv("data_in/phenotypes/COAT_WBC.csv")

CSFcols <- c(5:22, 24, 25)
WBCcols <- c(5:9)
invitrocols <- c(5:11)
###########################
#Read in variants files
###########################
classes<- c("character", "integer", rep("logical", 56),  "character", "character")
some93TF.snps <- read.csv("data_out/variants/df-some93_snpsREFvALT.csv", colClasses=classes)
some93otherTF.snps <- read.csv("data_out/variants/df-some93other_snpsREFvALT.csv", colClasses=classes)
some93TF.ind <- read.csv("data_out/variants/df-some93_indelsREFvALT.csv", colClasses=classes)
some93TF.ind$status <- "some93"
some93otherTF.ind <- read.csv("data_out/variants/df-some93other_indelsREFvALT.csv", colClasses=classes)
some93otherTF.ind$status <- "some93+other"

columnTypes <- c("character", "integer", "character", "character", "numeric", "integer", "integer", "character", "factor", "factor", "factor", "character", rep("character", 56), "character", "character")
some93.snps <- read.csv("data_out/variants/df-some93_snps.csv", colClasses=columnTypes)
some93other.snps <- read.csv("data_out/variants/df-some93other_snps.csv", colClasses=columnTypes) #685
some93.ind <- read.csv("data_out/variants/df-some93_indels.csv", colClasses=columnTypes)
some93other.ind <- read.csv("data_out/variants/df-some93other_indels.csv", colClasses=columnTypes) #685

some93TF.all <- rbind(some93TF.snps, some93otherTF.snps, some93TF.ind, some93otherTF.ind) #5605
some93.all <- rbind(some93.snps, some93other.snps, some93.ind, some93other.ind) #5605

####################
#Filter variants
####################
table(some93.all$HOM.VAR)
some93.4 <- subset(some93.all, HOM.VAR > 4 & HOM.VAR < 37) #1410 variants

inCEN <- c()
for(i in 1:length(some93.4[,1])){
  if(some93.4[i,]$POS > as.numeric(subset(centelo, V3 == "centromere" & V1 == paste0("chr", some93.4[i,]$CHROM))[1,]$V4) & some93.4[i,]$POS < as.numeric(subset(centelo, V1== paste0("chr", some93.4[i,]$CHROM))[1,]$V5)){
  inCEN <- append(inCEN, TRUE)
  }
  else inCEN <- append(inCEN, FALSE)
}

inTLO <- c()
for(i in 1:length(some93.4[,1])){
	tlo <- subset(centelo, V1== paste0("chr", some93.4[i,]$CHROM) & V3 == "telomere")
  if(some93.4[i,]$POS < as.numeric(tlo[1,]$V5) |  some93.4[i,]$POS > as.numeric(tlo[2,]$V4))  inTLO[i] <- TRUE
	else inTLO[i] <-  FALSE
}

some93.4$inCEN <- inCEN
some93.4$inTLO <- inTLO

some93.4NCT <- subset(some93.4, inCEN==FALSE & inTLO ==FALSE)
table(some93.4NCT$effect)
some93.4NCT$numALT <- stri_count_fixed(some93.4NCT$ALT, ",") #remove all variants with more than one alternative base (all had p = 1)
some93.4NCT <- subset(some93.4NCT, numALT ==0)

some93.eff <- subset(some93.4NCT, effect %nin% c("SYNONYMOUS_CODING", "INTERGENIC")) #781  variants, half of these are upstream
some93.eff  <- some93.eff[, c(1:12, 29:49, 51:53, 55:73)] #keep only COAT ST93
names(some93.eff)[13] <- "UgCl212"

#########################################
#Need to divide into CSF and WBC because
#missing some of the CSF data:
#missing  461 and 541
#########################################
names(some93.eff) %in% paste0("UgCl", dCSF$line)
some93.eff.CSF <- some93.eff[, names(some93.eff) %in% paste0("UgCl", dCSF$line)]
some93.eff.CSF <- cbind(some93.eff.CSF, some93.eff[,c(1:12, 51)])
some93.eff.CSF <- some93.eff.CSF[order(some93.eff.CSF$CP),]

nonVar.dCSF <- data.frame(some93.eff.CSF[1,])
var.dCSF <- data.frame(some93.eff.CSF[1,])
ks <- c()
for(k in 1:length(some93.eff.CSF[,1])){
	temp <- some93.eff.CSF[k, 1:35]
	alleleVarTable <- table(t(temp))
	alleleVarMin <- min(alleleVarTable)
	alleleVarMax <- max(alleleVarTable)
  if(alleleVarMax > c(sum(alleleVarTable)-4) | length(alleleVarTable)==1){
			ks <- append(ks, k)
			nonVar.dCSF <- rbind(nonVar.dCSF, some93.eff.CSF[k,])
		}
		else var.dCSF <- rbind(var.dCSF,  some93.eff.CSF[k,])
}
nonVar.dCSF <- nonVar.dCSF[-1,]
var.dCSF <- var.dCSF[-1,]

length(var.dCSF[,1]) #541 with  >4
write.csv(var.dCSF, "data_out/variants/var_dCSF.csv", row.names=FALSE)

########################
#CSF Stats
########################
j <- 0
estimCSF <- list()
pvalCSF.logR <- list()
odds.CSF <- list()
pstat.CSF <- list()
for(i in CSFcols){
	j <- j+1
	estims <- c()
	pvals.logR <- c()
	odds <- c()
	oddsRatio <- c()
	pStats <- c()
	for(k in 1:length(var.dCSF[,1])){
    pred <- t(dCSF[,i]) #phenotypic data
    resp <- t(var.dCSF[k,1:35]) #variant information
		dat <- data.frame(pred=t(pred), resp)
		names(dat) <- c("pred", "resp")
	  dat <- subset(dat, !is.na(pred))
	  dat <- subset(dat, !is.na(resp))
	  dat$resp <- as.character(dat$resp)
	  dat$resp[dat$resp!= var.dCSF[k, 38]] <- 1		#REF
	  dat$resp[dat$resp== var.dCSF[k, 38]] <- 0		#ALT
	  dat$pred <- as.numeric(dat$pred)
	  dat$resp <- as.numeric(dat$resp)
		t2 <- glm(resp~pred, data=dat, family=binomial)
		estims[k] <- summary(t2)$coefficients[2,1]
		pvals.logR[k] <- summary(t2)$coefficients[2,4]
		oddsRatio[k] <- logistic.display(t2)$table[1]
		pStats[k] <- paste0("est = ", round(summary(t2)$coefficients[2,1], 2), ", error = ", round(summary(t2)$coefficients[2,2], 2), ", z = ", round(summary(t2)$coefficients[2,3], 2), ", p = ", round(summary(t2)$coefficients[2,4], 3))
	}
	estimCSF[[j]] <- estims
	pvalCSF.logR[[j]] <- pvals.logR
	odds.CSF[[j]] <- oddsRatio
	pstat.CSF[[j]] <- pStats
	}
	names(estimCSF) <- names(dCSF)[c(13:29, 32)]
	names(pvalCSF.logR) <- names(dCSF)[c(13:29, 32)]
	names(odds.CSF) <- names(dCSF)[c(13:29, 32)]
	names(pstat.CSF) <- names(dCSF)[c(13:29, 32)]

nREF <- unlist(apply(var.dCSF, 1, function(x) table(x[1:35]==x["REF"])[2]))
nALT <- unlist(apply(var.dCSF, 1, function(x) table(x[1:35]==x["REF"])[1]))

geneP.CSF <- cbind(var.dCSF[, c(36, 37, 43, 44)], nREF, nALT, "IL1b" = unlist(pvalCSF.logR[[1]]), "IL2" = unlist(pvalCSF.logR[[2]]), "IL4" = unlist(pvalCSF.logR[[3]]), "IL5" = unlist(pvalCSF.logR[[4]]), "IL6" = unlist(pvalCSF.logR[[5]]), "IL7" = unlist(pvalCSF.logR[[6]]), "IL8" = unlist(pvalCSF.logR[[7]]), "IL10" = unlist(pvalCSF.logR[[8]]), "IL12" = unlist(pvalCSF.logR[[9]]), "IL13" = unlist(pvalCSF.logR[[10]]), "IL17" = unlist(pvalCSF.logR[[11]]), "GMCSF" = unlist(pvalCSF.logR[[12]]), "GCSF" = unlist(pvalCSF.logR[[13]]) , "IFNg" = unlist(pvalCSF.logR[[14]]), "MCP1" = unlist(pvalCSF.logR[[15]]), "MIP1b" = unlist(pvalCSF.logR[[16]]), "TNFa" = unlist(pvalCSF.logR[[17]]),"LFA_Titer" = unlist(pvalCSF.logR[[18]]), "protein" =  unlist(pvalCSF.logR[[19]]), var.dCSF[,c(38:42, 45:47)])
geneP.CSF$CP <- paste(geneP.CSF$CHROM, geneP.CSF$POS, sep=".")

geneOdds.CSF <- cbind(var.dCSF[, c(36, 37, 43, 44)], nREF, nALT, "IL1b" = unlist(odds.CSF[[1]]), "IL2" = unlist(odds.CSF[[2]]), "IL4" = unlist(odds.CSF[[3]]), "IL5" = unlist(odds.CSF[[4]]), "IL6" = unlist(odds.CSF[[5]]), "IL7" = unlist(odds.CSF[[6]]), "IL8" = unlist(odds.CSF[[7]]), "IL10" = unlist(odds.CSF[[8]]), "IL12" = unlist(odds.CSF[[9]]), "IL13" = unlist(odds.CSF[[10]]), "IL17" = unlist(odds.CSF[[11]]), "GMCSF" = unlist(odds.CSF[[12]]), "GCSF" = unlist(odds.CSF[[13]]) , "IFNg" = unlist(odds.CSF[[14]]), "MCP1" = unlist(odds.CSF[[15]]), "MIP1b" = unlist(odds.CSF[[16]]), "TNFa" = unlist(odds.CSF[[17]]),"LFA_Titer" = unlist(odds.CSF[[18]]), "protein" =  unlist(odds.CSF[[19]]), var.dCSF[,c(38:42, 45:47)])
geneOdds.CSF$CP <- paste(geneOdds.CSF$CHROM, geneOdds.CSF$POS, sep=".")

pstats.CSF <- cbind(var.dCSF[, c(36, 37, 43, 44)], nREF, nALT, "IL1b" = unlist(pstat.CSF[[1]]), "IL2" = unlist(pstat.CSF[[2]]), "IL4" = unlist(pstat.CSF[[3]]), "IL5" = unlist(pstat.CSF[[4]]), "IL6" = unlist(pstat.CSF[[5]]), "IL7" = unlist(pstat.CSF[[6]]), "IL8" = unlist(pstat.CSF[[7]]), "IL10" = unlist(pstat.CSF[[8]]), "IL12" = unlist(pstat.CSF[[9]]), "IL13" = unlist(pstat.CSF[[10]]), "IL17" = unlist(pstat.CSF[[11]]), "GMCSF" = unlist(pstat.CSF[[12]]), "GCSF" = unlist(pstat.CSF[[13]]) , "IFNg" = unlist(pstat.CSF[[14]]), "MCP1" = unlist(pstat.CSF[[15]]), "MIP1b" = unlist(pstat.CSF[[16]]), "TNFa" = unlist(pstat.CSF[[17]]),"LFA_Titer" = unlist(pstat.CSF[[18]]), "protein" =  unlist(pstat.CSF[[19]]), var.dCSF[,c(38:42, 45:47)])
pstats.CSF$CP <- paste(pstats.CSF$CHROM, pstats.CSF$POS, sep=".")

traitSig <- apply(geneP.CSF[,7:25], 2, function(x) length(subset(x, x < 0.05)))

geneP.CSF$numSig <- apply(geneP.CSF[,7:24], 1, function(x) length(subset(c(x), c(x) < 0.05)))
names(geneP.CSF)[3:4] <- c("gene", "effect")
names(geneP.CSF)[30:32] <- c("impact", "class", "AAchange")
write.csv(geneP.CSF, "data_out/GWAS/geneP-CSF.csv", row.names=FALSE)
write.csv(geneOdds.CSF, "data_out/GWAS/geneOdds-CSF.csv", row.names=FALSE)
write.csv(pstats.CSF, "data_out/GWAS/geneStats-CSF.csv", row.names=FALSE)

##########################
#WBC
##########################
some93.eff.WBC <- some93.eff[, names(some93.eff) %in% paste0("UgCl", dWBC$line)]
some93.eff.WBC <- cbind(some93.eff.WBC, some93.eff[,c(1:12, 51)])
some93.eff.WBC <- some93.eff.WBC[order(some93.eff.WBC$CP),]

nonVar.dWBC <- data.frame(some93.eff.WBC[1,])
var.dWBC <- data.frame(some93.eff.WBC[1,])
ks <- c()
for(k in 2:length(some93.eff.WBC[,1])){
	temp <- some93.eff.WBC[k, 15:47]
	alleleVarTable <- table(t(temp))
	alleleVarMin <- min(alleleVarTable)
	alleleVarMax <- max(alleleVarTable)
  if(alleleVarMax > c(sum(alleleVarTable)-4) | length(alleleVarTable)==1){
			ks <- append(ks, k)
			nonVar.dWBC <- rbind(nonVar.dWBC, some93.eff.WBC[k,])
		}
		else var.dWBC <- rbind(var.dWBC,  some93.eff.WBC[k,])
}
nonVar.dWBC <- nonVar.dWBC[-1,]
length(var.dWBC[,1]) #652 variants

write.csv(var.dWBC, "data_out/variants/var_dWBC.csv", row.names=FALSE)

j <- 0
estimWBC <- list()
pvalWBC.logR <- list()
odds.WBC <- list()
pstat.WBC <- list()
for(i in WBCcols){
	j <- j+1
	estims <- c()
	pvals.logR <- c()
	oddsRatio <- c()
	pStats <- c()
	for(k in 1:length(var.dWBC[,1])){
        pred <- t(dWBC[,i]) #phenotypic data
        resp <- t(var.dWBC[k,1:37]) #variant information
		dat <- data.frame(pred=t(pred), resp)
		names(dat) <- c("pred", "resp")
	    dat <- subset(dat, !is.na(pred))
	    dat <- subset(dat, !is.na(resp))
	    dat$resp <- as.character(dat$resp)
	    dat$resp[dat$resp!= var.dWBC[k, 40]] <- 1		#REF
	    dat$resp[dat$resp== var.dWBC[k, 40]] <- 0		#ALT
	    dat$pred <- as.numeric(dat$pred)
	    dat$resp <- as.numeric(dat$resp)
		t2 <- glm(resp~pred, data=dat, family=binomial)
		estims[k] <- summary(t2)$coefficients[2,1]
		pvals.logR[k] <- summary(t2)$coefficients[2,4]
	#	odds[k] <-  exp(coef(t2)[2])
		if(summary(t2)$coefficients[2,4] != 1) oddsRatio[k] <- logistic.display(t2)$table[1]
		else oddsRatio[k] <- NA
		pStats[k] <- paste0("est = ", round(summary(t2)$coefficients[2,1], 2), ", error = ", round(summary(t2)$coefficients[2,2], 2), ", z = ", round(summary(t2)$coefficients[2,3], 2), ", p = ", round(summary(t2)$coefficients[2,4], 3))
	    	}
	estimWBC[[j]] <- estims
	pvalWBC.logR[[j]] <- pvals.logR
	odds.WBC[[j]] <- oddsRatio
	pstat.WBC[[j]] <- pStats
	}
	names(estimWBC) <- names(dWBC)[WBCcols]
	names(pvalWBC.logR) <- names(dWBC)[WBCcols]
	names(odds.WBC) <- names(dWBC)[WBCcols]
	names(pstat.WBC) <- names(dWBC)[WBCcols]

 par(mfrow=c(1, 5), mar=c(1, 1, 1, 1))
 for(i in 1:5){
    hist( pvalWBC.logR[[i]], main=names(pvalWBC.logR)[i], xaxt="n", yaxt="n")
    abline(v=0.05, col="red")
}

nREF.WBC <- unlist(apply(var.dWBC, 1, function(x) table(x[1:37]==x["REF"])[2]))
nALT.WBC <- unlist(apply(var.dWBC, 1, function(x) table(x[1:37]==x["REF"])[1]))

geneP.WBC <- data.frame(var.dWBC[,c(38, 39, 45, 46)], nREF.WBC, nALT.WBC, "Survival" = unlist(pvalWBC.logR[[1]]), "hivrna" = unlist(pvalWBC.logR[[2]]), "cd4" = unlist(pvalWBC.logR[[3]]), "csf_wbc" = unlist(pvalWBC.logR[[4]]), "efa" = unlist(pvalWBC.logR[[5]]), var.dWBC[,c(40:44, 47:50)])
names(geneP.WBC)[3:4] <- c("gene", "effect")
names(geneP.WBC)[17:19] <- c("impact", "class", "AAchange")
geneP.WBC$CP <- paste(geneP.WBC$CHROM, geneP.WBC$POS, sep=".")

geneOdds.WBC <- data.frame(var.dWBC[,c(38, 39, 45, 46)], nREF.WBC, nALT.WBC, "Survival" = unlist(odds.WBC[[1]]), "hivrna" = unlist(odds.WBC[[2]]), "cd4" = unlist(odds.WBC[[3]]), "csf_wbc" = unlist(odds.WBC[[4]]), "efa" = unlist(odds.WBC[[5]]), var.dWBC[,c(40:44, 47:50)])
names(geneOdds.WBC)[3:4] <- c("gene", "effect")
names(geneOdds.WBC)[17:19] <- c("impact", "class", "AAchange")
geneOdds.WBC$CP <- paste(geneOdds.WBC$CHROM, geneOdds.WBC$POS, sep=".")

pstats.WBC <- data.frame(var.dWBC[,c(38, 39, 45, 46)], nREF.WBC, nALT.WBC, "Survival" = unlist(pstat.WBC[[1]]), "hivrna" = unlist(pstat.WBC[[2]]), "cd4" = unlist(pstat.WBC[[3]]), "csf_wbc" = unlist(pstat.WBC[[4]]), "efa" = unlist(pstat.WBC[[5]]), var.dWBC[,c(40:44, 47:50)])
names(pstats.WBC)[3:4] <- c("gene", "effect")
names(pstats.WBC)[17:19] <- c("impact", "class", "AAchange")
names(pstats.WBC)[17:19] <- c("impact", "class", "AAchange")
pstats.WBC$CP <- paste(pstats.WBC$CHROM, pstats.WBC$POS, sep=".")

traitSig.WBC <- apply(geneP.WBC[,7:11], 2, function(x) length(subset(x, x < 0.05)))
geneP.WBC$numSig <- apply(geneP.WBC[,7:11], 1, function(x) length(subset(c(x), c(x) < 0.05)))

write.csv(geneP.WBC, "data_out/GWAS/geneP-WBC.csv", row.names=FALSE)
write.csv(geneOdds.WBC, "data_out/GWAS/geneOdds-WBC.csv", row.names=FALSE)
write.csv(pstats.WBC, "data_out/GWAS/geneStats-WBC.csv", row.names=FALSE)


##########################
#invitro
##########################
some93.eff.invitro <- some93.eff[, names(some93.eff) %in% paste0("UgCl", strains$line[2:56])]
some93.eff.invitro <- cbind(some93.eff.invitro, some93.eff[,c(1:12, 51)])
some93.eff.invitro <- some93.eff.invitro[order(some93.eff.invitro$CP),]
var.dinvitro <-some93.eff.invitro
write.csv(var.dinvitro, "data_out/variants/var.dinvitro", row.names=FALSE)

nonVar.dINV <- data.frame(some93.eff.invitro[1,])
var.dINV <- data.frame(some93.eff.invitro[1,])
ks <- c()
for(k in 2:length(some93.eff.invitro[,1])){
	temp <- some93.eff.invitro[k, 15:47]
	alleleVarTable <- table(t(temp))
	alleleVarMin <- min(alleleVarTable)
	alleleVarMax <- max(alleleVarTable)
  if(alleleVarMax > c(sum(alleleVarTable)-4) | length(alleleVarTable)==1){
			ks <- append(ks, k)
			nonVar.dINV <- rbind(nonVar.invitro, some93.eff.invitro[k,])
		}
		else var.dINV <- rbind(var.dINV,  some93.eff.invitro[k,])
}
nonVar.dINV <- nonVar.dINV[-1,]
length(var.dINV[,1]) #652 variants

j <- 0
estimINV <- list()
pvalINV.logR <- list()
odds.INV <- list()
pstat.INV <- list()
for(i in invitrocols){
	j <- j+1
	estims <- c()
	pvals.logR <- c()
	oddsRatio <- c()
	pStats <- c()
	for(k in 1:length(var.dINV[,1])){
        pred <- t(strains[,i]) #phenotypic data
        resp <- t(var.dinvitro[k,1:38]) #variant information
		dat <- data.frame(pred=t(pred), resp)
		names(dat) <- c("pred", "resp")
	    dat <- subset(dat, !is.na(pred))
	    dat <- subset(dat, !is.na(resp))
	    dat$resp <- as.character(dat$resp)
	    dat$resp[dat$resp!= var.dinvitro[k, 40]] <- 1		#REF
	    dat$resp[dat$resp== var.dinvitro[k, 40]] <- 0		#ALT
	    dat$pred <- as.numeric(dat$pred)
	    dat$resp <- as.numeric(dat$resp)
		t2 <- glm(resp~pred, data=dat, family=binomial)
		estims[k] <- summary(t2)$coefficients[2,1]
		pvals.logR[k] <- summary(t2)$coefficients[2,4]
		if(summary(t2)$coefficients[2,3] > 0.0001) oddsRatio[k] <- logistic.display(t2)$table[1]
		if(summary(t2)$coefficients[2,3] < 0.0001) oddsRatio[k] <- NA
		pStats[k] <- paste0("est = ", round(summary(t2)$coefficients[2,1], 2), ", error = ", round(summary(t2)$coefficients[2,2], 2), ", z = ", round(summary(t2)$coefficients[2,3], 2), ", p = ", round(summary(t2)$coefficients[2,4], 3))
		#odds[k] <-  exp(coef(t2)[2])
	    	}
	estimINV[[j]] <- estims
	pvalINV.logR[[j]] <- pvals.logR
	odds.INV[[j]] <- oddsRatio
	pstat.INV[[j]] <- pStats
	}
	names(estimINV) <- names(strains)[c(8:14)]
	names(pvalINV.logR) <- names(strains)[c(8:14)]
	names(odds.INV) <- names(strains)[c(8:14)]
	names(pstat.INV) <- names(strains)[c(8:14)]

nREF.INV <- unlist(apply(var.dINV, 1, function(x) table(x[1:38]==x["REF"])[2]))
nALT.INV <- unlist(apply(var.dINV, 1, function(x) table(x[1:38]==x["REF"])[1]))

geneP.INV <- data.frame(var.dINV[,c(38, 39, 44, 45)], nREF.INV, nALT.INV, "uptake" = unlist(pvalINV.logR[[1]]), "adherence" = unlist(pvalINV.logR[[2]]), "chitin" = unlist(pvalINV.logR[[3]]), "growth" = unlist(pvalINV.logR[[4]]), "FLC" = unlist(pvalINV.logR[[5]]), "AMP" = unlist(pvalINV.logR[[6]]), "SERT" = unlist(pvalINV.logR[[7]]), var.dINV[,c(40:44, 47:50)])
names(geneP.INV)[3:4] <- c("gene", "effect")
names(geneP.INV)[17:19] <- c("impact", "class", "AAchange")
geneP.INV$CP <- paste(geneP.INV$CHROM, geneP.INV$POS, sep=".")

geneOdds.INV <- data.frame(var.dINV[,c(38, 39, 45, 46)], nREF.INV, nALT.INV, "uptake" = unlist(odds.INV[[1]]), "adherence" = unlist(odds.INV[[2]]), "chitin" = unlist(odds.INV[[3]]), "growth" = unlist(odds.INV[[4]]), "FLC" = unlist(odds.INV[[5]]), "AMP" = unlist(odds.INV[[6]]), "SERT" = unlist(odds.INV[[7]]), var.dINV[,c(40:44, 47:50)])
names(geneOdds.INV)[3:4] <- c("gene", "effect")
names(geneOdds.INV)[17:19] <- c("impact", "class", "AAchange")
geneOdds.INV$CP <- paste(geneOdds.INV$CHROM, geneOdds.INV$POS, sep=".")

geneStats.INV <- data.frame(var.dINV[,c(38, 39, 45, 46)], nREF.INV, nALT.INV, "uptake" = unlist(pstat.INV[[1]]), "adherence" = unlist(pstat.INV[[2]]), "chitin" = unlist(pstat.INV[[3]]), "growth" = unlist(pstat.INV[[4]]), "FLC" = unlist(pstat.INV[[5]]), "AMP" = unlist(pstat.INV[[6]]), "SERT" = unlist(pstat.INV[[7]]), var.dINV[,c(40:44, 47:50)])
names(geneStats.INV)[3:4] <- c("gene", "effect")
names(geneStats.INV)[17:19] <- c("impact", "class", "AAchange")
geneStats.INV$CP <- paste(geneStats.INV$CHROM, geneStats.INV$POS, sep=".")

traitSig.INV <- apply(geneP.INV[,7:13], 2, function(x) length(subset(x, x < 0.05)))

geneP.INV$numSig <- apply(geneP.INV[,7:13], 1, function(x) length(subset(c(x), c(x) < 0.05)))
write.csv(geneP.INV, "data_out/GWAS/geneP-INV.csv", row.names=FALSE)
write.csv(geneOdds.INV, "data_out/GWAS/geneOdds-INV.csv", row.names=FALSE)
write.csv(geneStats.INV, "data_out/GWAS/geneStats-INV.csv", row.names=FALSE)
