library(Hmisc)
library(here)
##########################
#Load strain and gene information
##########################
strains <- read.csv("data_in/general/180727strainInfo.csv")

genes <- read.csv("data_in/general/Cryptococcus_neoformans_H99_genes_FDB.bed", header=TRUE)
genes$chrom <- paste("chr", genes$chrom, sep="")
genes$length <- genes$chromEnd - genes$chromStart

##########################
#Read SNP tables
##########################
classes<- c("character", "integer", rep("logical", 56),  "character", "character")
all93.TF <- read.csv("tables_intermediate/variants/df-allST93_snpsREFvALT.csv", colClasses=classes)
all93other.TF <- read.csv("tables_intermediate/variants/df-allST93other_snpsREFvALT.csv", colClasses=classes)
some93.TF <- read.csv("tables_intermediate/variants/df-some93_snpsREFvALT.csv", colClasses=classes)
some93other.TF <- read.csv("tables_intermediate/variants/df-some93other_snpsREFvALT.csv", colClasses=classes)
otherST.TF <- read.csv( "tables_intermediate/variants/df-some93other_snpsREFvALT.csv", colClasses=classes)

columnTypes <- c("character", "integer", "character", "character", "numeric", "integer", "integer", "character", "factor", "factor", "factor", "character", rep("character", 56), "character", "character")
all93 <- read.csv("tables_intermediate/variants/df-allST93_snps.csv", colClasses=columnTypes)
all93other <- read.csv("tables_intermediate/variants/df-allST93other_snps.csv", colClasses=columnTypes)
some93 <- read.csv("tables_intermediate/variants/df-some93_snps.csv", colClasses=columnTypes)
some93other <- read.csv("tables_intermediate/variants/df-some93other_snps.csv", colClasses=columnTypes)
otherST <- read.csv( "tables_intermediate/variants/df-some93other_snps.csv", colClasses=columnTypes)

##########################
#Read IND tables
##########################
all93.TF.IND <- read.csv("tables_intermediate/variants/df-allST93_indelsREFvALT.csv", colClasses=classes)
all93other.TF.IND <- read.csv("tables_intermediate/variants/df-allST93other_indelsREFvALT.csv", colClasses=classes)
some93.TF.IND <- read.csv("tables_intermediate/variants/df-some93_indelsREFvALT.csv", colClasses=classes)
some93other.TF.IND <- read.csv("tables_intermediate/variants/df-some93other_indelsREFvALT.csv", colClasses=classes)
otherST.TF.IND <- read.csv( "tables_intermediate/variants/df-some93other_indelsREFvALT.csv", colClasses=classes)

all93.IND <- read.csv("tables_intermediate/variants/df-allST93_indels.csv", colClasses=columnTypes)
all93other.IND <- read.csv("tables_intermediate/variants/df-allST93other_indels.csv", colClasses=columnTypes)
some93.IND <- read.csv("tables_intermediate/variants/df-some93_indels.csv", colClasses=columnTypes)
some93other.IND <- read.csv("tables_intermediate/variants/df-some93other_indels.csv", colClasses=columnTypes)
otherST.IND <- read.csv( "tables_intermediate/variants/df-some93other_indels.csv", colClasses=columnTypes)


#######################################
#look at number of variants x gene length
#######################################
all93.all <- rbind(all93, all93.IND)
all93.gene <- subset(all93, effect %in% c("NON_SYNONYMOUS_CODING", "STOP_GAINED", "START_GAINED", "CODON_CHANGE_PLUS_CODON_DELETION", "CODON_CHANGE_PLUS_CODON_INSERTION", "CODON_DELETION", "CODON_INSERTION", "FRAME_SHIFT", "SPLICE_SITE_ACCEPTOR", "SPLICE_SITE_DONOR", " STOP_LOST")) #2714 variants
all93.numVar.table <- table(all93.gene$gene)
all93.numVar <- data.frame(all93.numVar.table)
names(all93.numVar)[1] <- "name"

genes.all93 <- subset(genes, name %in% all93.gene$gene)
genes.all93 <- genes.all93[order(genes.all93$name),]
all93.numVar <- all93.numVar[order(all93.numVar$name),]
all93.numVar$chrom <- genes.all93$chrom
all93.numVar$chromStart <- genes.all93$chromStart
all93.numVar$chromEnd <- genes.all93$chromEnd
all93.numVar$length <- genes.all93$length
all93.numVar$alias <- genes.all93$alias
all93.numVar$description <- genes.all93$description

all93.syn <- subset(all93, effect == "SYNONYMOUS_CODING")
all93.numVar.table.syn <- table(all93.syn$gene)
all93.numVar.syn <- data.frame(all93.numVar.table.syn)
names(all93.numVar.syn)[1] <- "name"

genes.all93.syn <- subset(genes, name %in% all93.syn$gene)
genes.all93.syn <- genes.all93.syn[order(genes.all93.syn$name),]
all93.numVar.syn <- all93.numVar.syn[order(all93.numVar.syn$name),]
all93.numVar.syn$chrom <- genes.all93.syn$chrom
all93.numVar.syn$chromStart <- genes.all93.syn$chromStart
all93.numVar.syn$chromEnd <- genes.all93.syn$chromEnd
all93.numVar.syn$length <- genes.all93.syn$length
all93.numVar.syn$alias <- genes.all93.syn$alias
all93.numVar.syn$description <- genes.all93.syn$description

all93.numVar.syn.unique <- subset(all93.numVar.syn, name %nin% all93.numVar$name) #909 genes have synonymous variants with no genic variants
all93.numVar.unique <- subset(all93.numVar, name %nin% all93.numVar.syn$name) #879 genes have genic variants without synonymous variants
extra.all93.syn <- data.frame(all93.numVar.syn.unique$name, Freq=rep(0, length(all93.numVar.syn.unique$name)), chrom=all93.numVar.syn.unique$chrom, chromStart=all93.numVar.syn.unique$chromStart, chromEnd = all93.numVar.syn.unique$chromEnd, length=all93.numVar.syn.unique$length, all93.numVar.syn.unique$alias, all93.numVar.syn.unique$description)
names(extra.all93.syn) <- names(all93.numVar.syn)
all93.numVar.genic <- rbind(all93.numVar, extra.all93.syn)
all93.numVar.genic$chrom <- factor(all93.numVar.genic$chrom, levels=c(paste0("chr", 1:14)))
all93.numVar.genic <- all93.numVar.syn[order(all93.numVar.genic$chrom, all93.numVar.genic$chromStart),]

extra.all93.genic <- data.frame(all93.numVar.unique$name, Freq=rep(0, length(all93.numVar.unique$name)), chrom=all93.numVar.unique$chrom, chromStart=all93.numVar.unique$chromStart, chromEnd = all93.numVar.unique$chromEnd, length=all93.numVar.unique$length, all93.numVar.unique$alias, all93.numVar.unique$description)
names(extra.all93.genic) <- names(all93.numVar.syn)
all93.numVar.syn <- rbind(all93.numVar.syn, extra.all93.genic)
all93.numVar.syn$chrom <- factor(all93.numVar.syn$chrom, levels=c(paste0("chr", 1:14)))
all93.numVar.syn <- all93.numVar.syn[order(all93.numVar.syn$chrom, all93.numVar.syn$chromStart),]

numGenesChr.ddn <- split(all93.numVar.syn, all93.numVar.syn$chrom)
numGenesChr <- unlist(lapply(numGenesChr.ddn, function(x) length(x[,1])))
numGenesChrCS <- cumsum(numGenesChr)

plot(seq_along(all93.numVar.genic[,1]), all93.numVar.genic$Freq, col="blue", xaxt="n")
points(seq_along(all93.numVar.syn[,1]), all93.numVar.syn$Freq, col="red")
axis(1, at= numGenesChrCS, labels=FALSE)
text(c(numGenesChr[1]/2, numGenesChr[2]/2+numGenesChrCS[1], numGenesChr[3]/2+numGenesChrCS[2],  numGenesChr[4]/2+numGenesChrCS[3], numGenesChr[5]/2+numGenesChrCS[4], numGenesChr[6]/2+numGenesChrCS[5], numGenesChr[7]/2+numGenesChrCS[6], numGenesChr[8]/2+numGenesChrCS[7], numGenesChr[9]/2+numGenesChrCS[8], numGenesChr[10]/2+numGenesChrCS[9], numGenesChr[11]/2+numGenesChrCS[10], numGenesChr[12]/2+numGenesChrCS[11], numGenesChr[13]/2+numGenesChrCS[12], numGenesChr[14]/2+numGenesChrCS[13]), -10, c(paste0("chr", 1:14)), xpd=NA, cex=0.7)
