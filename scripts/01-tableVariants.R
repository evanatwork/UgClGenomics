
library(Hmisc)
##########################
#Load strain information
##########################
strains <- read.csv("data_in/general/180727strainInfo.csv")

genes <- read.csv("data_in/general/Cryptococcus_neoformans_H99_genes_FDB.bed", header=TRUE, as.is=TRUE)
genes$chrom <- paste("chr", genes$chrom, sep="")
genes$length <- genes$chromEnd - genes$chromStart
genes.CNAG <- genes[grep("CNAG", genes$name),]

################################################################################
#SNPS
################################################################################
dfall <- read.csv("data_in/gatk_processed/UgClSeq_snps.txt", as.is=TRUE)
#For each potential variant count how many of the genomes do not have any base information that positions
table(dfall$numDat)
dfall <- subset(dfall, numDat >49) #127344
names(dfall)[8:12] <- c("gene", "effect", "impact", "class", "aaChange")

#Simplify effects
dfall$effect[dfall$effect=="SYNONYMOUS_STOP"] <- "SYNONYMOUS_CODING"
dfall$effect[dfall$effect=="SPLICE_SITE_ACCEPTOR"] <- "SPLICE_SITE"
dfall$effect[dfall$effect=="SPLICE_SITE_DONOR"] <- "SPLICE_SITE"
dfall$effect <- factor(dfall$effect, levels=c("SYNONYMOUS_CODING", "NON_SYNONYMOUS_CODING", "START_GAINED","STOP_GAINED", "UPSTREAM", "DOWNSTREAM", "UTR_5_PRIME", "UTR_3_PRIME", "NON_SYNONYMOUS_START", "INTERGENIC", "INTRON", NA))

#Score REF vs. ALT counts
snps.all <- data.frame(CP=dfall$CP, numDat = dfall$numDat)
for(i in 13:68){
  pos <- as.character(dfall[,i])
  #this means that NA will sort of be scored as alternative, which is maybe not ideal. We can't score on ALT because sometimes there are two options.
  is.ref <- pos == as.character(dfall$REF)
  snps.all <- cbind(snps.all, is.ref)
  }
names(snps.all)[3:58] <- names(dfall)[13:68]

#calculate the number of non-reference genomes for each variant
nonREF <- unlist(apply(snps.all[,3:58], 1, function(x) table(t(x))["FALSE"]))
snps.all$nonREF <- nonREF
snps.all$nonREF[is.na(snps.all$nonREF)] <- 0

#Subset into different dataframes
allALT <- subset(snps.all, nonREF==56) #2941
names(allALT) <- names(snps.all)
allALT$status <- "allAlt"
variable <- subset(snps.all, nonREF!= 0 & nonREF!=56) #124403
names(variable) <- names(snps.all)

ST93 <- c(3, 4, 16, 18:39, 41:43, 45:58)
notST93 <- c(5:15, 17, 40, 44)
status <- c()
for(i in 1:length(variable[,1])) {
  if(FALSE %in% variable[i, ST93]){
    if(TRUE %nin% variable[i, ST93]){
      if(FALSE %nin% variable[i, notST93]) status[i] <- "all93"
      if(FALSE %in% variable[i, notST93]) status[i] <- "all93+other"
    }
    if(TRUE %in% variable[i, ST93]){
      if(FALSE %nin% variable[i, notST93]) status[i] <- "some93"
      if(FALSE %in% variable[i, notST93]) status[i] <- "some93+other"
    }
   }
 if(FALSE %nin% variable[i, ST93]){
   if(FALSE %in% variable[i, notST93]) status[i] <- "otherST"
 }
}
variable$status <- status

names(variable)
#Continue subset
all93 <- subset(variable, status=="all93") #4681
all93other <- subset(variable, status=="all93+other") #40695
some93 <- subset(variable, status=="some93") #3396
some93other <- subset(variable, status=="some93+other") #1277
otherST <- subset(variable, status=="otherST") #74354

#write.csv(all93, "data_out/variants/df-allST93_snpsREFvALT.csv", row.names=FALSE)
#write.csv(all93other, "data_out/variants/df-allST93other_snpsREFvALT.csv", row.names=FALSE)
#write.csv(some93, "data_out/variants/df-some93_snpsREFvALT.csv", row.names=FALSE)
#write.csv(some93other, "data_out/variants/df-some93other_snpsREFvALT.csv", row.names=FALSE)
#write.csv(otherST, "data_out/variants/df-otherST_snpsREFvALT.csv", row.names=FALSE)

length(unique(dfall$gene)) #7509
all93.snps <- subset(dfall, CP %in% all93$CP)
length(unique(all93.snps$gene)) #2573
all93other.snps <- subset(dfall, CP %in% all93other$CP)
length(unique(all93other.snps$gene)) #5949
some93.snps <- subset(dfall, CP %in% some93$CP)
length(unique(some93.snps$gene)) #2185
some93other.snps <- subset(dfall, CP %in% some93other$CP)
length(unique(some93other.snps$gene)) #228
otherST.snps <- subset(dfall, CP %in% otherST$CP)
length(unique(otherST.snps$gene)) #7246

#write.csv(all93.snps, "data_out/variants/df-allST93_snps.csv", row.names=FALSE)
#write.csv(all93other.snps, "data_out/variants/df-allST93other_snps.csv", row.names=FALSE)
#write.csv(some93.snps, "data_out/variants/df-some93_snps.csv", row.names=FALSE)
#write.csv(some93other.snps, "data_out/variants/df-some93other_snps.csv", row.names=FALSE)
#write.csv(otherST.snps, "data_out/variants/df-otherST_snps.csv", row.names=FALSE)

################################################################################
#INDELS
################################################################################
dfall.IND <- read.csv("data_in/gatk_processed/UgClSeq_indels.txt")

#For each potential variant count how many of the genomes do not have any base information that positions
dfall.IND <- subset(dfall.IND, numDat >49) #15032
names(dfall.IND)[8:12] <- c("gene", "effect", "impact", "class", "aaChange")

#Score REF vs. ALT counts
indels.all <- data.frame(CP=dfall.IND$CP, numDat = dfall.IND$numDat)
for(i in 13:68){
  pos <- as.character(dfall.IND[,i])
  #this means that NA will sort of be scored as alternative, which is maybe not ideal. We can't score on ALT because sometimes there are two options.
  is.ref <- pos == as.character(dfall.IND$REF)
  indels.all <- cbind(indels.all, is.ref)
  }
names(indels.all)[3:58] <- names(dfall.IND)[13:68]

#calculate the number of non-reference genomes for each variant
nonREF.IND <- unlist(apply(indels.all[,3:58], 1, function(x) table(t(x))["FALSE"]))
indels.all$nonREF <- nonREF.IND
indels.all$nonREF[is.na(indels.all$nonREF)] <- 0

#Subset into different dataframes
allALT.IND <- subset(indels.all, nonREF==56) #301
names(allALT.IND) <- names(indels.all)
allALT.IND$status <- "allAlt"
variable.IND <- subset(indels.all, nonREF!= 0 & nonREF!=50) #14704
names(variable.IND) <- names(indels.all)

status.IND <- c()
for(i in 1:length(variable.IND[,1])) {
  if(FALSE %in% variable.IND[i, ST93]){
    if(TRUE %nin% variable.IND[i, ST93]){
      if(FALSE %nin% variable.IND[i, notST93]) status.IND[i] <- "all93"
      if(FALSE %in% variable.IND[i, notST93]) status.IND[i] <- "all93+other"
    }
    if(TRUE %in% variable.IND[i, ST93]){
      if(FALSE %nin% variable.IND[i, notST93]) status.IND[i] <- "some93"
      if(FALSE %in% variable.IND[i, notST93]) status.IND[i] <- "some93+other"
    }
   }
 if(FALSE %nin% variable.IND[i, ST93]){
   if(FALSE %in% variable.IND[i, notST93]) status.IND[i] <- "otherST"
 }
}
variable.IND$status.IND <- status.IND
names(variable.IND)[60] <- "status"

#Continue subset
all93.IND <- subset(variable.IND, status=="all93") #429
all93other.IND <- subset(variable.IND, status=="all93+other") #4838
some93.IND <- subset(variable.IND, status=="some93") #479
some93other.IND <- subset(variable.IND, status=="some93+other") #453
otherST.IND <- subset(variable.IND, status=="otherST") #8505

#write.csv(all93.IND, "data_out/variants/df-allST93_indelsREFvALT.csv", row.names=FALSE)
#write.csv(all93other.IND, "data_out/variants/df-allST93other_indelsREFvALT.csv", row.names=FALSE)
#write.csv(some93.IND, "data_out/variants/df-some93_indelsREFvALT.csv", row.names=FALSE)
#write.csv(some93other.IND, "data_out/variants/df-some93other_indelsREFvALT.csv", row.names=FALSE)
#write.csv(otherST.IND, "data_out/variants/df-otherST_indelsREFvALT.csv", row.names=FALSE)

all93.indels <- subset(dfall.IND, CP %in% all93.IND$CP) #12483
all93other.indels <- subset(dfall.IND, CP %in% all93other.IND$CP) #13904
some93.indels <- subset(dfall.IND, CP %in% some93.IND$CP) #4122
some93other.indels <- subset(dfall.IND, CP %in% some93other.IND$CP) #685
otherST.indels <- subset(dfall.IND, CP %in% otherST.IND$CP) #29970

#write.csv(all93.indels, "data_out/variants/df-allST93_indels.csv", row.names=FALSE)
#write.csv(all93other.indels, "data_out/variants/df-allST93other_indels.csv", row.names=FALSE)
#write.csv(some93.indels, "data_out/variants/df-some93_indels.csv", row.names=FALSE)
#write.csv(some93other.indels, "data_out/variants/df-some93other_indels.csv", row.names=FALSE)
#write.csv(otherST.indels, "data_out/variants/df-otherST_indels.csv", row.names=FALSE)

######################
#all variants
######################
allVar <- rbind(dfall, dfall.IND)
#write.csv(allVar, "data_out/variants/df-allVar.csv", row.names=FALSE)
allVar <- read.csv("data_out/variants/df-allVar.csv", as.is=TRUE)

nrow(allVar)
length(unique(allVar$gene)) #7561 genes with at least one mutations

sub_genes <- subset(genes, name %in% allVar$gene)
nrow(subset(sub_genes, description == "hypothetical RNA")) #690
nrow(subset(sub_genes, description == "hypothetical protein")) #2982
tRNAs <- sub_genes[grep("-tRNA", sub_genes$description),] #75
#with descriptions = 7561 - 690 - 2982 = 3889
#Percent without description = 2982/(2982+3889)


table(allVar$effect)

#upstream
nrow(subset(allVar, effect %in% c("UPSTREAM", "UTR_5_PRIME")))/nrow(allVar) # 0.428
#downstream
nrow(subset(allVar, effect %in% c("DOWNSTREAM", "UTR_3_PRIME")))/nrow(allVar) #0.106
#synonymous
nrow(subset(allVar, effect %in% c("SYNONYMOUS_CODING")))/nrow(allVar) #0.221
nrow(subset(allVar, effect %in% c("INTERGENIC")))/nrow(allVar) #.028


allVar.gene <- subset(allVar, effect %in% c("NON_SYNONYMOUS_CODING", "STOP_GAINED", "START_GAINED", "CODON_CHANGE_PLUS_CODON_DELETION", "CODON_CHANGE_PLUS_CODON_INSERTION", "CODON_DELETION", "CODON_INSERTION", "FRAME_SHIFT", "SPLICE_SITE_ACCEPTOR", "SPLICE_SITE_DONOR", " STOP_LOST"))
nrow(allVar.gene)/nrow(allVar)
length(unique(allVar.gene$gene))
nrow(subset(allVar.gene, effect %in% c("NON_SYNONYMOUS_CODING")))/nrow(allVar.gene) #0.896

allVar.table <- table(allVar$gene)
allVar.gene.table <- table(allVar.gene$gene)

#add in the genes that don't have any variants
noVar <- genes$name %nin% names(allVar.table)
noVar <- genes$name[noVar] #945 genes have no identified variants - but some of these are tRNAs
allGenes.table <- c(allVar.table, rep(0, length(noVar)))
names(allGenes.table) <- c(names(allVar.table), noVar)
allGenes.table <- allGenes.table[order(names(allGenes.table))]
genes <- genes[order(genes$name),]
genes$numVar <- allGenes.table
genes.CNAG <- genes[grep("CNAG", genes$name),]
t <- table(genes.CNAG$description)
gene.descript  <- subset(genes, description %nin% c("hypothetical protein", "hypothetical RNA")) #4256 genes with a description



############################
#FIGURES
############################
pdf("manuscript/figures/Figure1B-180727numVar-genLength.pdf", width=4, height=4)
plot(gene.descript$length, gene.descript$numVar, xlab="gene length", ylab="number of variants", yaxt="n")
axis(2, las=2)
abline(lm(gene.descript$numVar~gene.descript$length), col="red")
dev.off()

cor.test(gene.descript$numVar, gene.descript$length)
#t = 33.001, df = 4254, p-value < 2.2e-16
#0.4514712

noVar.CNAG <- noVar[grep("CNAG", noVar)] #778 genes have no identified variants

allGenes.table <- c(allVar.table, rep(0, length(noVar)))
names(allGenes.table) <- c(names(allVar.table), noVar)

allGenes.table <- allGenes.table[order(names(allGenes.table))]
genes <- genes[order(genes$name),]
genes$numVar <- allGenes.table

pdf("manuscript/figures/Figure1A-180727numVar-all.pdf", width=6, height=4)
par(fig=c(0, 1, 0, 1))
hist(genes.CNAG$numVar, xaxt="n", yaxt="n", ylab="number of genes", xlab="number of variants", breaks=100, ylim=c(0, 1400), main="")
axis(1, pos=0)
axis(2, las=2, pos=0)
par(fig = c(0.3,0.9, 0.3, 0.95), new = T)
hist(subset(genes.CNAG, numVar > 50)$numVar, yaxt="n", xaxt="n", ylab="", xlab="", breaks=100, main="")
axis(1, pos=0)
axis(2, las=2, pos=50)
dev.off()

length(subset(genes.CNAG, numVar > 50)$numVar) #435

#How many variants per genome?
snps <- c()
for (i in 13:68){
	sub <- subset(dfall, dfall[,i] != ".")
	snps[i-12] <- table(as.character(sub[,i]) == as.character(sub$REF))[1]
}

indels <- c()
for (i in 13:68){
	sub <- subset(dfall.IND, dfall.IND[,i] != ".")
	indels[i-12] <- table(as.character(sub[,i]) == as.character(sub$REF))[1]
}

#characterize by strain
strains$numSNP <- snps
strains$numIND <- indels
strains$numVAR <- indels+snps
ord.numVAR <- strains[order(strains$numVAR),]
ord.numVAR <- ord.numVAR[c(1:12, 22, 13:21, 23:56),]

pdf("manuscript/figures/Figure1C-180727numSNPs.pdf", width=9.5, height=5)
par(oma=c(1, 1, 1, 4))
t <- barplot(t(cbind(ord.numVAR$numIND, ord.numVAR$numSNP)), yaxt="n", ylim=c(0, 60000))
axis(2, las=2, pos=-1, cex.axis=0.8)
text(t-0.9, -2000, ord.numVAR$line, xpd=NA, srt=-45, cex=0.6, pos=4)
#text(t-0.9, -2000, paste0("ST", ord.numVAR$ST, "-", ord.numVAR$line), xpd=NA, srt=-45, cex=0.6, pos=4)
#mtext("Sequence Type - Strain name", side=1, line=2)
mtext("UgCl strain number", side=1, line=2)
mtext("Number of variants per strain", side=2, line=3)
par(xpd=NA)
 legend(50,50000,c("SNPs", "indels"), pch = c(22,22), col=c(grey(0.7), grey(0.3)), pt.bg=c(grey(0.7), grey(0.3)), cex=0.8, bg="white")
dev.off()
