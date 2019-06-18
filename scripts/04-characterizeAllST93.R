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

hypoRNAs <- read.csv("data_in/general/190422hypotheticalRNAs.csv")
locs <- substr(hypoRNAs$Genomic.Location, 11, 30)
locs.ddn <- strsplit(locs, " ")
chromStart <- as.numeric(unlist(lapply(locs.ddn, function(x) x[[1]])))
chromEnd <- as.numeric(unlist(lapply(locs.ddn, function(x) x[[3]])))
chrom <- paste0("chr", hypoRNAs$Chromosome)
Strand <- ifelse(hypoRNAs$Gene.Strand == "forward", "+", "-")
name <- hypoRNAs$Gene.ID
length <- chromEnd - chromStart
description <- "hypothetical RNA"
alias <- "null"
hypoRNAs.df <- data.frame(chrom, chromStart, chromEnd, length, Strand, name, description, alias)

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
all93.unnamed_genes <- subset(all93, gene %nin% genes$name) #357

all93.named_genes$type <- "temp"
all93.named_genes$type[all93.named_genes$effect == "SYNONYMOUS_CODING"] <- "syn"
all93.named_genes$type[all93.named_genes$effect %in% c("UTR_5_PRIME", "UTR_3_PRIME", "UPSTREAM", "DOWNSTREAM")] <-
"upstream"
all93.named_genes$type[all93.named_genes$type == "temp"] <- "ingene"

all93.named.ddn <- split(all93.named_genes, all93.named_genes$gene) #2573

names <- c()
freq.gene <- c()
freq.up <- c()
freq.syn <- c()
freq.other <- c()
total <-c()
for (i in 1:length(all93.named.ddn)){
  names <- append(names, as.character(names(all93.named.ddn)[i]))
  temp <- all93.named.ddn[[i]]
  freq.gene <- append(freq.gene, nrow(subset(temp, type == "ingene")))
  freq.syn <- append(freq.syn, nrow(subset(temp, type == "syn")))
  freq.up <- append(freq.up, nrow(subset(temp, type == "upstream")))
  }


all93.named.df <- data.frame(name= names, freq.gene, freq.up, freq.syn) #141

all93.named.df$chrom <- subset(genes, name %in% all93.named.df$name)$chrom
all93.named.df$chromStart <- subset(genes, name %in% all93.named.df$name)$chromStart
all93.named.df$chromEnd <- subset(genes, name %in% all93.named.df$name)$chromEnd
all93.named.df$length <- subset(genes, name %in% all93.named.df$name)$length
all93.named.df$alias <- subset(genes, name %in% all93.named.df$name)$alias
all93.named.df$description <- subset(genes, name %in% all93.named.df$name)$description
all93.named.df$Freq <- all93.named.df$freq.up + all93.named.df$freq.gene + all93.named.df$freq.syn

##############################################################
#Add frequencies of different types into master dataframe
##############################################################
all93_wFreq <- all93.named.df[, c("name", "chrom", "chromStart", "chromEnd", "length", "alias", "description", "freq.gene", "freq.up", "freq.syn", freq.gene = "Freq")]
all93_wFreq$chrom <- factor(all93_wFreq$chrom, levels=paste0("chr", 1:14))


all93_wFreq <- all93_wFreq[order(all93_wFreq$chrom, all93_wFreq$chromStart),] #put back into position order
all93_wFreq_high <-  subset(all93_wFreq, all93_wFreq$Freq > 10)
write.csv(all93_wFreq_high, "manuscript/tables/TableS2_all93-highvariants.csv", row.names=FALSE)

numGenesChr.ddn <- split(all93_wFreq, all93_wFreq$chrom)
numGenesChr <- unlist(lapply(numGenesChr.ddn, function(x) length(x[,1])))
numGenesChrCS <- cumsum(numGenesChr)


#pdf("manuscript/figures/Figure3_ST93-numVar-gene.pdf", width=12, height=4)
tiff(filename = "manuscript/figures/Figure3_ST93-numVar-gene.tiff", width = 5, height = 2.5, units = 'in', res = 300, compression = 'lzw', pointsize = 9)
par(mar=c(1,1,1,1), oma=c(3, 4, 1, 1), fig=c(0, 1, 0, 1), mgp=c(1,0.75,0))
plot(seq_along(all93_wFreq[,1]), all93_wFreq$Freq, col="black", xaxt="n", ylim=c(0, 50), yaxt="n", xlab="", ylab="Number of variants")
axis.break(axis=2,breakpos=45,pos=NA,bgcol="white",breakcol="black",  style="slash",brw=0.02)
points(22, 50)
axis(2, las=2, at=c(0, 10, 20, 30, 40, 49), labels=c(0, 10, 20, 30, 40,  90))
axis(1, at= numGenesChrCS, labels=FALSE)
text(c(numGenesChr[1]/2, numGenesChr[2]/2+numGenesChrCS[1], numGenesChr[3]/2+numGenesChrCS[2],  numGenesChr[4]/2+numGenesChrCS[3], numGenesChr[5]/2+numGenesChrCS[4], numGenesChr[6]/2+numGenesChrCS[5], numGenesChr[7]/2+numGenesChrCS[6], numGenesChr[8]/2+numGenesChrCS[7], numGenesChr[9]/2+numGenesChrCS[8], numGenesChr[10]/2+numGenesChrCS[9], numGenesChr[11]/2+numGenesChrCS[10], numGenesChr[12]/2+numGenesChrCS[11], numGenesChr[13]/2+numGenesChrCS[12], numGenesChr[14]/2+numGenesChrCS[13]), -5, c(1:14), xpd=NA, cex=0.7)
mtext("Position in genome (chromosome)", side=1, line=2)
mtext("Number of variants", side=2, line=2.5)
dev.off()

########################################################################
#Need unnamed genes too - these are the hypothetical RNAs and intergenic
########################################################################
all93.unnamed_genes <- subset(all93, gene %nin% genes$name) #357
all93.isna <- subset(all93.unnamed_genes, is.na(gene)) #138
all93.rna <- subset(all93.unnamed_genes, !is.na(gene)) #219

all93.rna.ddn <- split(all93.rna, all93.rna$gene) #141

names <- c()
freq.gene <- c()
freq.up <- c()
freq.syn <- c()
freq.other <- c()
total <-c()
for (i in 1:length(all93.rna.ddn)){
  names <- append(names, as.character(names(all93.unnamed_genes.ddn)[i]))
  temp <- all93.rna.ddn[[i]]
  freq.up <- append(freq.up, nrow(subset(temp, effect %in% c("UTR_5_PRIME", "UTR_3_PRIME", "UPSTREAM", "DOWNSTREAM"))))
  freq.syn <- append(freq.syn, nrow(subset(temp, effect %in% c("SYNONYMOUS_CODING"))))
  freq.gene <- append(freq.gene, nrow(subset(temp, effect %nin% c("UTR_5_PRIME", "UTR_3_PRIME", "UPSTREAM", "DOWNSTREAM","SYNONYMOUS_CODING"))))
}

all93.rna.df <- data.frame(name= names, freq.gene, freq.up, freq.syn) #141
all93.rna.df.tRNA <- subset(all93.rna.df, names == "CNAG_10023")
all93.rna.df <- subset(all93.rna.df, name !="CNAG_10023")
all93.rna.df$chrom <- subset(hypoRNAs.df, name %in% all93.rna.df$name)$chrom
all93.rna.df$chromStart <- subset(hypoRNAs.df, name %in% all93.rna.df$name)$chromStart
all93.rna.df$chromEnd <- subset(hypoRNAs.df, name %in% all93.rna.df$name)$chromEnd
all93.rna.df$length <- subset(hypoRNAs.df, name %in% all93.rna.df$name)$length
all93.rna.df$alias <- subset(hypoRNAs.df, name %in% all93.rna.df$name)$alias
all93.rna.df$description <- subset(hypoRNAs.df, name %in% all93.rna.df$name)$description
all93.rna.df$Freq <- all93.rna.df$freq.up + all93.rna.df$freq.gene + all93.rna.df$freq.syn #140 rows
#218

all93_gene_RNAs <- rbind(all93_wFreq, all93.rna.df) #2713
all93_gene_RNAs$name <- as.character(all93_gene_RNAs$name)
all93_gene_RNAs$chrom <- as.character(all93_gene_RNAs$chrom)
all93_gene_RNAs$alias <- as.character(all93_gene_RNAs$alias)
all93_gene_RNAs$description <- as.character(all93_gene_RNAs$description)

all93.rna.df.tRNA <- c(name= "CNAG_10023", chrom = "chr2", chromStart = "911318", chromEnd = "911427", length = 109, alias= "NULL", description = "tRNA-OTHER",  freq.gene = 0, freq.up = 0, freq.syn = 0, Freq = 1)

all93_gene_RNAs <- rbind(all93_gene_RNAs, all93.rna.df.tRNA)
all93_gene_RNAs$Freq <- as.numeric(all93_gene_RNAs$Freq)

all93.intergenic <- data.frame(all93.isna$CHROM)
names(all93.intergenic)[1] <- "chrom"
all93.intergenic$name <- NA
all93.intergenic$chromStart <- all93.isna$POS
all93.intergenic$chromEnd <- NA
all93.intergenic$length <- 1
all93.intergenic$alias <- NA
all93.intergenic$description <- "intergenic"
all93.intergenic$freq.gene <- 0
all93.intergenic$freq.up <- 0
all93.intergenic$freq.syn <- 0
all93.intergenic$Freq <- 1

all93_gene_RNAs_inter <- rbind(all93_gene_RNAs, all93.intergenic)
all93_gene_RNAs_inter <- all93_gene_RNAs_inter[order(all93_gene_RNAs_inter$name),]
names(all93_gene_RNAs_inter)[11] <- "total"
write.csv(all93_gene_RNAs_inter, "manuscript/tables/TableS1-all93variants.csv", row.names=FALSE)
