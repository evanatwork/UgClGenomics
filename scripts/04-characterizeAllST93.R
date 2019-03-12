library(Hmisc)
library(here)
library(plotrix)

##########################
#Load strain information
##########################
strains <- read.csv("data_in/general/180727strainInfo.csv")

###########################
#Read in files
###########################
allSeq.snps <- read.csv("data_in/gatk_processed/UgClSeq_snps.txt", as.is=TRUE)
allSeq.indels <- read.csv("data_in/gatk_processed/UgClSeq_indels.txt")
allSeq.snps$type <- "snps"
allSeq.indels$type <- "indels"

genes <- read.csv("data_in/general/180709genesList.csv")
genes$length <- genes$chromEnd - genes$chromStart
genes <- genes[order(genes$name),]
genes.CNAG <- genes[grep("CNAG", genes$name),]

chrLens <- c(2291499, 1621675, 1575141, 1084805, 1814975, 1422463, 1399503, 1398693, 1186808, 1059964, 1561994, 774062, 756744, 926563)

centelo <- read.table("data_in/general/centromeres_telomeres_cna3.gff", sep="\t", as.is=TRUE)
centelo <- rbind(centelo, c("chr5", "prediction", "MAT", 185700 , 275730, ".", "+", ".", "locus_tag Chr5 MAT"))

###########################
#Characterize the whole set
###########################
allSeq <- rbind(allSeq.snps, allSeq.indels)
allSeq <- subset(allSeq, numDat > 40) #144542
names(allSeq)[8:12] <- c("gene", "effect", "impact", "class", "aaChange")
table(allSeq$type)

numVar.UgCl <- c()
for(i in 13:52){
  numVar.UgCl[i-12] <- table(allSeq[,i]==allSeq$REF)["FALSE"]
}

#what type of variant?
table(allSeq$effect)

#remove "SYNONYMOUS_CODING", "SYNONYMOUS_STOP", "INTERGENIC"
allSeq <- subset(allSeq, effect %nin% c("SYNONYMOUS_CODING", "SYNONYMOUS_STOP", "INTERGENIC"))

length(unique(allSeq$gene)) #7500

######################################################
#Characterize the set that's in all of the ST93 genomes
######################################################
#read in data
columnTypes <- c("character", "integer", "character", "character", "numeric", "integer", "integer", "character", "factor", "factor", "factor", "character", rep("character", 56), "character", "character")

all93.snps <- read.csv("tables_intermediate/variants/df-allST93_snps.csv", colClasses=columnTypes) #4681
all93.indels <- read.csv("tables_intermediate/variants/df-allST93_indels.csv", colClasses=columnTypes) #429

#join SNPs and INDELS
table(all93.snps$gene %in% genes$name) #4353 true, 328 false
table(all93.indels$gene %in% genes$name) #400 true, 29 false
all93 <- rbind(all93.snps, all93.indels)  #5110
all93.named_genes <- subset(all93, gene %in% genes$name) #4753

all93.numVar <- table(all93$gene)
all93.df <- data.frame(all93.numVar)
names(all93.df)[1] <- "name"
all93.df$name <- as.character(all93.df$name) #2714
all93.df.named_genes <- subset(all93.df, name %in% genes$name)
#write.csv(all93.df, "manuscript/tables/TableS1-allST93variants.csv", row.names=FALSE)

# #Look at the variants in genes (not hypothetical RNAs)
# genes.all93 <- subset(genes, name %in% all93.df$name)
# genes.all93 <- genes.all93[order(genes.all93$name),]
#
# all93.gene.df <- subset(all93.df, name %in% genes$name) #2573
# all93.gene.df <- all93.gene.df[order(all93.gene.df$name),]
# all93.gene.df$chrom <- genes.all93$chrom
# all93.gene.df$chrom <- factor(all93.gene.df$chrom, levels=c(paste0("chr", 1:14)))
# all93.gene.df$chromStart <- genes.all93$chromStart
# all93.gene.df$chromEnd <- genes.all93$chromEnd
# all93.gene.df$length <- genes.all93$length
# all93.gene.df$alias <- genes.all93$alias
# all93.gene.df$description <- genes.all93$description
#
#
# #190309- THIS IS NOT IN THE MS ATM
# all93.not_gene.df <- subset(all93.df, name %nin% genes$name) #141

######################################################
#create dataframe that counts the number of allST93 variants in genes
######################################################
all93.gene <- subset(all93.named_genes, effect %in% c("NON_SYNONYMOUS_CODING", "STOP_GAINED", "START_GAINED", "CODON_CHANGE_PLUS_CODON_DELETION", "CODON_CHANGE_PLUS_CODON_INSERTION", "CODON_DELETION", "CODON_INSERTION", "FRAME_SHIFT", "SPLICE_SITE_ACCEPTOR", "SPLICE_SITE_DONOR", " STOP_LOST")) #1267
all93.gene.numVar <- table(all93.gene$gene)
all93.gene.df <- data.frame(all93.gene.numVar) #994/988
names(all93.gene.df)[1] <- "name"
all93.gene.df$name <- as.character(all93.gene.df$name)

genes.all93 <- subset(genes, name %in% all93.gene.df$name)
genes.all93 <- genes.all93[order(genes.all93$name),]
all93.gene.df <- all93.gene.df[order(all93.gene.df$name),]
genes.all93$Freq <- all93.gene.df$Freq

all93.gene.df_extra <- subset(all93.df.named_genes, name %nin% all93.gene.df$name)  #1585
nogenes.all93 <- subset(genes, name %in% all93.gene.df_extra$name)
nogenes.all93$Freq <- 0

all93.gene.full <- rbind(genes.all93, nogenes.all93)
all93.gene.full <- all93.gene.full[order(all93.gene.full$name),]
##############################################################
#create dataframe that counts the number of allST93 synonymous variants
##############################################################
all93.syn <- subset(all93.named_genes, effect == "SYNONYMOUS_CODING")
all93.syn.numVar <- table(all93.syn$gene)
all93.syn.df <- data.frame(all93.syn.numVar)
names(all93.syn.df)[1] <- "name"
all93.syn.df$name <- as.character(all93.syn.df$name)

genes.all93.syn <- subset(genes, name %in% all93.syn.df$name)
genes.all93.syn <- genes.all93.syn[order(genes.all93.syn$name),]
all93.syn.df <- all93.syn.df[order(all93.syn.df$name),] #794
genes.all93.syn$Freq <- all93.syn.df$Freq

all93.syn.df_extra <- subset(all93.df.named_genes, name %nin% all93.syn.df$name)  #1779
nogenes.all93.syn <- subset(genes, name %in% all93.syn.df_extra$name)
nogenes.all93.syn$Freq <- 0

all93.syn.full <- rbind(genes.all93.syn, nogenes.all93.syn)
all93.syn.full <- all93.syn.full[order(all93.syn.full$name),]


######################################################
#create dataframe that counts the number of allST93 variants upstream or downstream - but nearly all upstream
######################################################
all93.up <- subset(all93.named_genes, effect %in% c("UTR_5_PRIME", "UTR_3_PRIME", "UPSTREAM", "DOWNSTREAM")) #2443 variants
all93.up.numVar <- table(all93.up$gene) #1469
all93.up.df <- data.frame(all93.up.numVar)
names(all93.up.df)[1] <- "name"
all93.up.df$name <- as.character(all93.up.df$name)

genes.all93.up <- subset(genes, name %in% all93.up.df$name)
genes.all93.up <- genes.all93.up[order(genes.all93.up$name),]
all93.up.df <- all93.up.df[order(all93.up.df$name),]
genes.all93.up$Freq <- all93.up.df$Freq

all93.up.df_extra <- subset(all93.df.named_genes, name %nin% all93.up.df$name)  #1104
nogenes.all93.up <- subset(genes, name %in% all93.up.df_extra$name)
nogenes.all93.up$Freq <- 0

all93.up.full <- rbind(genes.all93.up, nogenes.all93.up)
all93.up.full <- all93.up.full[order(all93.up.full$name),]
##############################################################
#Add frequencies of different types into master dataframe
##############################################################
all93_wFreq <- all93.gene.full[, c("name", "chrom", "chromStart", "chromEnd", "length", "alias", "description", "freq.gene" = "Freq")]
names(all93_wFreq)[8] <- "freq.gene"
all93_wFreq$freq.up <- all93.up.full$Freq
all93_wFreq$freq.syn <- all93.syn.full$Freq
all93_wFreq <- all93_wFreq[order(all93_wFreq$chrom, all93_wFreq$chromStart),] #put back into position order
#FISX
all93_wFreq$Freq <- all93_wFreq$freq.gene + all93_wFreq$freq.up + all93_wFreq$freq.syn

write.csv(all93_wFreq, "manuscript/tables/TableS1-all93variants.csv", row.names=FALSE)

all93_wFreq_high <-  subset(all93_wFreq, all93_wFreq$Freq > 10)
write.csv(all93_wFreq_high, "manuscript/tables/TableS2_all93-highvariants.csv", row.names=FALSE)


numGenesChr.ddn <- split(all93_wFreq, all93_wFreq$chrom)
numGenesChr <- unlist(lapply(numGenesChr.ddn, function(x) length(x[,1])))
numGenesChrCS <- cumsum(numGenesChr)


pdf("manuscript/figures/Figure3_ST93-numVar-gene.pdf", width=12, height=4)
plot(seq_along(all93.df[,1]), all93.df$Freq, col="black", xaxt="n", ylim=c(0, 50), yaxt="n", xlab="", ylab="Number of variants")
axis.break(axis=2,breakpos=45,pos=NA,bgcol="white",breakcol="black",  style="slash",brw=0.02)
points(22, 50)
axis(2, las=2, at=c(0, 10, 20, 30, 40, 49), labels=c(0, 10, 20, 30, 40,  90))
axis(1, at= numGenesChrCS, labels=FALSE)
text(c(numGenesChr[1]/2, numGenesChr[2]/2+numGenesChrCS[1], numGenesChr[3]/2+numGenesChrCS[2],  numGenesChr[4]/2+numGenesChrCS[3], numGenesChr[5]/2+numGenesChrCS[4], numGenesChr[6]/2+numGenesChrCS[5], numGenesChr[7]/2+numGenesChrCS[6], numGenesChr[8]/2+numGenesChrCS[7], numGenesChr[9]/2+numGenesChrCS[8], numGenesChr[10]/2+numGenesChrCS[9], numGenesChr[11]/2+numGenesChrCS[10], numGenesChr[12]/2+numGenesChrCS[11], numGenesChr[13]/2+numGenesChrCS[12], numGenesChr[14]/2+numGenesChrCS[13]), -5, c(paste0("chr", 1:14)), xpd=NA, cex=0.7)
mtext("Position in genome", side=1, line=2)
dev.off()


