library(Hmisc)
library(stringr)
library(epiDisplay)

lenU <- function(x) length(unique(x))

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

genes <- read.csv("data_in/general/180709genesList.csv", as.is=TRUE)

###########################
#Read in phenotype files
###########################
dCSF <- read.csv("data_in/phenotypes/Clean_Compiled_CSF_log2duplicates_d0.csv", as.is=TRUE) #"duplicates", t=0
dCSF$IL1b <- as.numeric(dCSF$IL1b)

dWBC <- read.csv("data_in/phenotypes/COAT_WBC.csv")

###########################
#Read in variants files
###########################
var.dCSF <- read.csv("tables_intermediate/variants/var_dCSF.csv")
var.dCSF$CP <- paste(var.dCSF$CHROM, var.dCSF$POS, sep=".")
var.dWBC <- read.csv("tables_intermediate/variants/var_dWBC.csv")
var.dWBC$CP <- paste(var.dWBC$CHROM, var.dWBC$POS, sep=".")

###########################
#Read in significant variants files
###########################
geneP.CSF <- read.csv("tables_intermediate/GWAS/geneP-CSF.csv")
geneP.WBC <- read.csv("tables_intermediate/GWAS/geneP-WBC.csv")
geneP.INV <- read.csv("tables_intermediate/GWAS/geneP-INV.csv")

###########################
#Read in logistic regression variants files
###########################
geneOdds.CSF <- read.csv("tables_intermediate/GWAS/geneOdds-CSF.csv", colClasses ="character")
geneOdds.CSF$CP <- paste(geneOdds.CSF$CHROM, geneOdds.CSF$POS, sep=".")
geneOdds.WBC <- read.csv("tables_intermediate/GWAS/geneOdds-WBC.csv", colClasses = "character")
geneOdds.WBC$CP <- paste(geneOdds.WBC$CHROM, geneOdds.WBC$POS, sep=".")
geneOdds.INV <- read.csv("tables_intermediate/GWAS/geneOdds-INV.csv", colClasses = "character")
geneOdds.INV$CP <- paste(geneOdds.INV$CHROM, geneOdds.INV$POS, sep=".")

###########################
#Read in stats info files
###########################
geneStats.CSF <- read.csv("tables_intermediate/GWAS/geneStats-CSF.csv", colClasses ="character")
geneStats.CSF$CP <- paste(geneStats.CSF$CHROM, geneStats.CSF$POS, sep=".")
geneStats.WBC <- read.csv("tables_intermediate/GWAS/geneStats-WBC.csv", colClasses = "character")
geneStats.WBC$CP <- paste(geneStats.WBC$CHROM, geneStats.WBC$POS, sep=".")
geneStats.INV <- read.csv("tables_intermediate/GWAS/geneStats-INV.csv", colClasses = "character")
geneStats.INV$CP <- paste(geneStats.INV$CHROM, geneStats.INV$POS, sep=".")

##########################
#Common variants
##########################
common.WBC<- subset(geneP.WBC, CP  %in% geneP.CSF$CP)
common.INV<- subset(geneP.INV, CP  %in% geneP.CSF$CP)
geneP.all <- cbind(geneP.CSF[,1:25], common.WBC[,7:11], common.INV[,7:13], geneP.CSF[,26:34])
numSig <- apply(geneP.all[,7:37], 1, function(x) length(subset(c(x), c(x) < 0.05)))
geneP.all$numSig <- numSig
geneP.all$effect <- as.character(geneP.all$effect)
geneP.all$effect[geneP.all$effect=="NON_SYNONYMOUS_CODING"] <- "ns"
geneP.all$effect[geneP.all$effect=="START_GAINED"] <- "start+"
geneP.all$effect[geneP.all$effect=="UTR_5_PRIME"] <- "UTR-5"
geneP.all$effect[geneP.all$effect=="UTR_3_PRIME"] <- "UTR-3"
geneP.all$effect[geneP.all$effect=="UPSTREAM"] <- "upstream"
geneP.all$effect[geneP.all$effect=="FRAME_SHIFT"] <- "frameshift"
geneP.all$effect[geneP.all$effect=="DOWNSTREAM"] <- "downstream"

alias <- c()
for(i in 1:length(geneP.all$gene)){
  if(is.na(as.character(geneP.all$gene[i]))) alias[i] <- "null"
	if (as.character(geneP.all$gene[i]) %nin% genes$name) alias[i]<- "null"
  else alias[i] <- as.character(subset(genes, name==as.character(geneP.all$gene[i]))$alias)
}
geneP.all$alias <- alias
write.csv(geneP.all, "tables_intermediate/GWAS/180808genePall.csv", row.names=FALSE)

geneP.sig.all <- subset(geneP.all, numSig > 0)
geneP.sig.2 <- subset(geneP.sig.all, numSig > 1)
multiVar <- names(subset(table(geneP.sig.all$gene), table(geneP.sig.all$gene)>1))
geneP.sig.multi <- subset(geneP.sig.all, gene %in% multiVar)
geneP.sig.2.multi <- subset(geneP.sig.2,  CP %in% geneP.sig.multi$CP)
geneP.sig.2.multi$class <- "ab"
geneP.sig.2.only <- subset(geneP.sig.2,  CP %nin% geneP.sig.multi$CP)
geneP.sig.2.only$class <- "a"
geneP.sig.multi.only <- subset(geneP.sig.multi,  CP %nin% geneP.sig.2$CP)
geneP.sig.multi.only$class <- "b"
geneP.sig.allTraits <- rbind(geneP.sig.multi.only, geneP.sig.2.only, geneP.sig.2.multi)
geneP.sig.allTraits$CP <- paste(geneP.sig.allTraits$CHROM, geneP.sig.allTraits $POS, sep=".")
geneP.sig.allTraits $gene <- as.character(geneP.sig.allTraits$gene)

#WBC CP traits only
WBC.P <- subset(geneP.WBC, CP %nin% geneP.CSF$CP)
inVitro.P <- subset(geneP.INV, CP %nin% geneP.CSF$CP)

WBC_inVitro.P <- cbind(WBC.P[,1:11], inVitro.P[,7:13], WBC.P[,12:20])
numSig <- apply(WBC_inVitro.P[,7:18], 1, function(x) length(subset(c(x), c(x) < 0.05)))
WBC_inVitro.P$numSig <- numSig

WBC_inVitro.P$effect <- as.character(WBC_inVitro.P$effect)
WBC_inVitro.P$effect[WBC_inVitro.P$effect=="NON_SYNONYMOUS_CODING"] <- "ns"
WBC_inVitro.P$effect[WBC_inVitro.P$effect=="START_GAINED"] <- "start+"
WBC_inVitro.P$effect[WBC_inVitro.P$effect=="UTR_5_PRIME"] <- "UTR-5"
WBC_inVitro.P$effect[WBC_inVitro.P$effect=="UTR_3_PRIME"] <- "UTR-3"
WBC_inVitro.P$effect[WBC_inVitro.P$effect=="UPSTREAM"] <- "upstream"
WBC_inVitro.P$effect[WBC_inVitro.P$effect=="FRAME_SHIFT"] <- "frameshift"
WBC_inVitro.P$effect[WBC_inVitro.P$effect=="DOWNSTREAM"] <- "downstream"

alias <- c()
for(i in 1:length(WBC_inVitro.P$gene)){
  if(is.na(as.character(WBC_inVitro.P$gene[i]))) alias[i] <- "null"
	if (as.character(WBC_inVitro.P$gene[i]) %nin% genes$name) alias[i]<- "null"
  else alias[i] <- as.character(subset(genes, name==as.character(WBC_inVitro.P$gene[i]))$alias)
}
WBC_inVitro.P$alias <- alias

WBC_inVitro.sig <- subset(WBC_inVitro.P, numSig > 0)

#Gene is already in the all traits set
WBC_inVitro.sig.allTraits <- subset(WBC_inVitro.sig, gene %in% geneP.sig.allTraits$gene)
WBC_inVitro.sig.allTraits.2 <-  subset(WBC_inVitro.sig.allTraits, numSig > 1)
WBC_inVitro.sig.allTraits.2$class <- c("ab")
WBC_inVitro.sig.allTraits.1 <-  subset(WBC_inVitro.sig.allTraits, numSig ==1)
WBC_inVitro.sig.allTraits.1$class <- "b"

#Gene is not in the all traits set
WBC_inVitro.sig.unique <- subset(WBC_inVitro.sig, gene %nin% geneP.sig.allTraits$gene)
WBC_inVitro.sig.unique.2 <- subset(WBC_inVitro.sig.unique, numSig > 1)
WBC_inVitro.sig.unique.2$class <- "a"
WBC_inVitro.sig.unique.1 <- subset(WBC_inVitro.sig.unique, numSig == 1)

WBC_inVitro.sig.allTraits <- rbind(WBC_inVitro.sig.allTraits.2, WBC_inVitro.sig.allTraits.1, WBC_inVitro.sig.unique.2)

#save with each phenotype from each variant on a different row
sub <- c()
gene <-c()
chr <- c()
pos <- c()
pos1 <- c()
effect <- c()
class <- c()
pheno <- c()
numSigSum <- c()
numVar <- c()
numPhen <- c()
oddsProb <- c()
stats <- c()
for(i in unique(geneP.sig.allTraits$gene)){
	sub <- subset(geneP.sig.allTraits, gene==i)
	allP <- unlist(sub[,7:37])
	numSig <- subset(allP, allP < 0.05)
		for(j in 1: nrow(sub)){
			temp <- sub[j, 7:37]
			temp.CP <- sub[j, "CP"]
			sigP <- names(sub)[7:37][which(temp < 0.05)]
			if(length(sigP)==1){
					gene <- append(gene, sub$gene[j])
					chr <- append(chr, sub$CHROM[j])
					pheno <- append(pheno, paste(sigP, collapse=";"))
					effect <- append(effect, sub[j, "effect"])
					class <- append(class, sub[j, "class"])
					pos <- append(pos, sub[j, "POS"])
					numPhen <- append(numPhen, length(sigP))
					if(sigP %in% names(geneOdds.CSF)) oddsProb <- append(oddsProb, subset(geneOdds.CSF, CP==temp.CP) [sigP])
					if(sigP %in% names(geneOdds.WBC)) oddsProb <- append(oddsProb, subset(geneOdds.WBC, CP==temp.CP) [sigP])
					if(sigP %in% names(geneOdds.INV)) oddsProb <- append(oddsProb, subset(geneOdds.INV, CP==temp.CP) [sigP])
					if(sigP %in% names(geneStats.CSF)) stats <- append(stats, subset(geneStats.CSF, CP==temp.CP) [sigP])
					if(sigP %in% names(geneStats.WBC)) stats <- append(stats, subset(geneStats.WBC, CP==temp.CP) [sigP])
					if(sigP %in% names(geneStats.INV)) stats <- append(stats, subset(geneStats.INV, CP==temp.CP) [sigP])
				}
			if(length(sigP)>1){
				 oddsTemp <- c()
				 statsTemp <- c()
					for(k in 1:length(sigP)){
						gene <- append(gene, sub$gene[j])
						chr <- append(chr, sub$CHROM[j])
						pheno <- append(pheno, sigP[k])
						effect <- append(effect, sub[j, "effect"])
						class <- append(class, sub[j, "class"])
						pos <- append(pos, sub[j, "POS"])
						if(sigP[k] %in% names(geneOdds.CSF)) oddsTemp[k] <-  subset(geneOdds.CSF, CP==temp.CP) [sigP[k]]
						if(sigP[k] %in% names(geneOdds.WBC))  oddsTemp[k] <-  subset(geneOdds.WBC, CP==temp.CP) [sigP[k]]
						if(sigP[k] %in% names(geneOdds.INV))  oddsTemp[k]  <- subset(geneOdds.INV, CP==temp.CP) [sigP[k]]
						if(sigP[k] %in% names(geneStats.CSF)) statsTemp[k] <- subset(geneStats.CSF, CP==temp.CP) [sigP[k]]
						if(sigP[k] %in% names(geneStats.WBC)) statsTemp[k] <- subset(geneStats.WBC, CP==temp.CP) [sigP[k]]
						if(sigP[k] %in% names(geneStats.INV)) statsTemp[k] <- subset(geneStats.INV, CP==temp.CP) [sigP[k]]
					}
				oddsProb <- append(oddsProb, oddsTemp)
				stats <- append(stats, statsTemp)
		}
	}
}

allSig_allPhen <- data.frame(gene, chr, pos, class, effect, pheno, oddsProb= unlist(oddsProb), stats= unlist(stats))
allSig_allPhen$CP <- paste(allSig_allPhen$chr, allSig_allPhen$pos, sep=".")
allSig_allPhen <- allSig_allPhen[order(allSig_allPhen$gene),]
allSig_allPhen <- allSig_allPhen[, c("gene", "chr", "pos", "class", "effect", "pheno", "oddsProb", "stats")]
names(allSig_allPhen)[7] <- "odds ratio (CI)"
names(allSig_allPhen)[8] <- "logistic regression"
#write.csv(allSig_allPhen, "manuscript/tables/TableS7_geneP2-allsig-stats.csv", row.names=FALSE)

#save with each variant on a different row
sub <- c()
gene <-c()
chr <- c()
pos <- c()
pos1 <- c()
effect <- c()
class <- c()
pheno <- c()
numSigSum <- c()
numVar <- c()
numPhen <- c()
for(i in unique(geneP.sig.allTraits$gene)){
	sub <- subset(geneP.sig.allTraits, gene==i)
	allP <- unlist(sub[,7:37])
	numSig <- subset(allP, allP < 0.05)
		for(j in 1: nrow(sub)){
			temp <- sub[j, 7:37]
			sigP <- names(sub)[7:37][which(temp < 0.05)]
			if(length(sigP)!=0){
					gene <- append(gene, sub$gene[j])
					chr <- append(chr, sub$CHROM[j])
					pheno <- append(pheno, paste(sigP, collapse=";"))
					effect <- append(effect, sub[j, "effect"])
					class <- append(class, sub[j, "class"])
					pos <- append(pos, sub[j, "POS"])
					numPhen <- append(numPhen, length(sigP))
		}
	}
}

allSig <- data.frame(gene, chr, pos, class, effect, pheno, numPhen)
allSig$gene <- as.character(allSig$gene)
allSig$pheno <- as.character(allSig$pheno)
allSig <- allSig[order(allSig$gene),] #133 positions in 41 genes

sub <- c()
gene <-c()
chr <- c()
pos <- c()
pos1 <- c()
effect <- c()
pheno <- c()
numSigSum <- c()
numVar <- c()
class <- c()
numPhen <- c()
for(i in unique(WBC_inVitro.sig.allTraits$gene)){
	sub <- subset(WBC_inVitro.sig.allTraits, gene==i)
	allP <- unlist(sub[, 7:18])
	for(j in 1: length(sub[,1])){
			temp <- sub[j, 7:18]
			sigP <- names(sub)[7:18][which(temp < 0.05)]
			if(length(sigP)!=0){
				# if(length(sigP) >1){
					gene <- append(gene, as.character(sub$gene[j]))
					chr <- append(chr, sub$CHROM[j])
					pheno <- append(pheno, as.character(paste(sigP, collapse=";")))
					effect <- append(effect, sub[j, "effect"])
					class <- append(class, sub[j, "class"])
					pos <- append(pos, sub[j, "POS"])
					numPhen <- append(numPhen, length(sigP))
				}
			}
		}
# 	}
# }

allSig.WBC_inVitro <- data.frame(gene, chr, pos,class, effect, pheno, numPhen)
allSig.WBC_inVitro$gene <- as.character(allSig.WBC_inVitro$gene)
allSig.WBC_inVitro$pheno <- as.character(allSig.WBC_inVitro$pheno)

allSigCP <- rbind(allSig, allSig.WBC_inVitro)
allSigCP <- allSigCP[order(allSigCP$gene),] #145 positions in 42 genes

alias <- c()
for(i in 1:length(allSigCP$gene)){
  if(is.na(as.character(allSigCP$gene[i]))) alias[i] <- "null"
  else alias[i] <- as.character(subset(genes, name==as.character(allSigCP$gene[i]))$alias)
}
allSigCP$alias <- alias
allSigCP$numSig <- allSigCP$numSigSum

subset(allSigCP, gene %nin% genes$name)
#CNAG_12610, CNAG_13108, CNAG_13204
#from fungidb, all are hypothetical RNAs. Add manually into genes file based on fungidb information
genes <- rbind(genes, c(6998, "chr7", 45069, 45328, "-", "CNAG_12610", 259, "hypothetical RNA"), c(6999,"chr13", 128964, 129469, "+", "CNAG_13108",  505, "hypothetical RNA"), c(7000,"chr14", 918503, 919365, "-", "CNAG_13204",  862, "hypothetical RNA"))
genes$chromStart <- as.numeric(genes$chromStart)
genes$chromEnd <- as.numeric(genes$chromEnd)
genes$length <- as.numeric(genes$length)

allSigCP$loc <- locCalc(allSigCP)

write.csv(allSigCP, "manuscript/tables/TableS6_geneP-allSigCP.csv", row.names=FALSE)

