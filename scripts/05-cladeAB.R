library(here)
library(Hmisc)
library(remotes)
if (system.file(package="flexplot") == "") {
  library(lme4)
  remotes::install_github("dustinfife/flexplot")
}
if (system.file(package="fifer") == "") {
  library(mice)
  remotes::install_github("dustinfife/fifer")
}
library(fifer)

# #Plot variant locations function
snpLocLine <- function(list, list2, col1 = "purple", col2 = "orange", col3 = "purple",  listOrient="single"){
	j <- 0
	k <- 0
	par(mfrow=c(14, 1), mar=c(0.5, 1, 1, 1), oma=c(4, 1, 1, 1))
	for(i in 1:14){
		plot(c(1, chrLens[i]), c(1, 1), type="l", xlim=c(0, 2301499), xaxt="n", yaxt="n")
		if(i %in% names(list)){
			j <- j+1
			if(length(list[[j]][,1]) > 0) {
				listTemp0 <- list[[j]]
				if(listOrient=="both") arrows(listTemp0[, "POS"], c(1.3, 0.7), listTemp0[, "POS"], c(1.05, 0.95), length=0.025, col=col1, lwd=1.25)
				else arrows(listTemp0[, "POS"], c(1.3), listTemp0[, "POS"], c(1.05), length=0.025, col=col1, lwd=1.25)
			}
		}
		if(!missing(list2)){
		  if(i %in% names(list2)){
				k <- k+1
				if(length(list2[[k]][,1]) > 0) {
					listTemp <- list2[[k]]
					if(listOrient=="both") arrows(listTemp[, "POS"], c(1.3, 0.7), listTemp[, "POS"], c(1.05, 0.95), length=0.025, col=col1, lwd=1.25)
					else arrows(listTemp[, "POS"], c(0.7), listTemp[, "POS"], c(0.95), length=0.025, col=col2, lwd=1.25)
		    }
			}
	}
	CENTLO <- subset(centelo, V1==paste0("chr", i))
		for(m in 1:length(CENTLO[,1])){
			if(CENTLO[m,3]=="centromere"){
				points(c(CENTLO[m, 4],CENTLO[m, 5]), c(1, 1), type="l", col="grey", lwd=4)
			}
			if(CENTLO[m,3]=="telomere"){
				points(c(CENTLO[m, 4],CENTLO[m, 5]), c(1, 1), type="l", col="grey", lwd=4)
			}
	#		if(centelos[k,3]=="MAT"){
	#			points(c(centelos[k, 4],centelos[k, 5]), c(1, 1), type="l", col="grey", lwd=5)
	#		}
		}
		if(i ==14) axis(1, cex=0.8, at =c(0, 500000, 1000000, 1500000, 2000000), labels=c(0, 500, 1000, 1500, 2000))
		else axis(1, labels=FALSE)
		mtext(paste("Chr", i), side=3, adj=0.01, cex=0.8)
	}
	mtext("Position on chromosome (kb)", side=1, outer=TRUE, line=2)
}


##########################
#Load strain information
##########################
strains <- read.csv("data_in/general/180727strainInfo.csv")
strains <- subset(strains, ST=="93" & arm != "pre-COAT")

centelo <- read.table("data_in/general/centromeres_telomeres_cna3.gff", sep="\t", as.is=TRUE)
centelo <- rbind(centelo, c("chr5", "prediction", "MAT", 185700 , 275730, ".", "+", ".", "locus_tag Chr5 MAT"))

chrLens <- c(2291499, 1621675, 1575141, 1084805, 1814975, 1422463, 1399503, 1398693, 1186808, 1059964, 1561994, 774062, 756744, 926563)

genes <- read.csv("data_in/general/Cryptococcus_neoformans_H99_genes_FDB.bed", header=TRUE)
genes$chrom <- paste("chr", genes$chrom, sep="")

###########################
#Read in files
###########################
all93.snps <- read.csv("data_out/variants/df-allST93_snps.csv") #4681
all93.indels <- read.csv("data_out/variants/df-allST93_indels.csv") #429
all93 <- rbind(all93.snps, all93.indels)  #5110

classesTF<- c("character", "integer", rep("logical", 56),  "character", "character")
columnTypes <- c("character", "integer", "character", "character", "numeric", "integer", "integer", "character", "factor", "factor", "factor", "character", rep("character", 56), "character", "character")
some93TF <- read.csv("data_out/variants/df-some93_snpsREFvALT.csv", colClasses=classesTF) #3396
some93.snps <- read.csv("data_out/variants/df-some93_snps.csv", colClasses=columnTypes)
some93TF.ind <- read.csv("data_out/variants/df-some93_indelsREFvALT.csv", colClasses=classesTF) #479
some93.ind <- read.csv("data_out/variants/df-some93_indels.csv", colClasses=columnTypes)

###########################
#Characterize the whole set
###########################
some93TF.all <- rbind(some93TF, some93TF.ind)
some93.all <- rbind(some93.snps, some93.ind) #3875

#######################################################
#ID the SNPs and INDELs that are associated with ST93 clade A/B
#######################################################
some93 <- subset(some93.all, HOM.VAR > 13 & HOM.VAR < 24) #only 252 variants
some93.sub <- some93[, c(1:12, 29:49, 51:53, 55:68)] #Only keep UgCl COAT ST93 strains
some93.sub$CP <- paste(some93.sub$CHROM, some93$POS, sep=".")

baseTable <- c(CP =some93.sub[1,3], ST93A.ALT = 0, ST93A.REF = 0, ST93B.ALT = 0, ST93B.REF = 0)
for(i in 1:length(some93.sub[,1])){
	t <- data.frame(strain = names(some93.sub[13:50]), base = t(some93.sub[i, 13:50] ), ST93 = strains$ST93)
	REF <- some93.sub[i,]$REF
	names(t)[2]  <- "base"
	t$REF <- t$base == REF
	baseTable <- rbind(baseTable, c(CP = as.character(some93.sub[i,"CP"]), ST93A.ALT = as.numeric(table(t$ST93, t$REF)[1,1]), ST93A.REF = as.numeric(table(t$ST93, t$REF)[1,2]), ST93B.ALT = as.numeric(table(t$ST93, t$REF)[2,1]), ST93B.REF = as.numeric(table(t$ST93, t$REF)[2,2])))
 }

baseTable <- baseTable[-1,]
baseTable <- as.data.frame(baseTable)
baseTable$chr <- some93$CHROM
baseTable$pos <- some93$POS
baseTable$gene <- some93$SNPEFF_GENE_NAME
baseTable$effect <- some93$SNPEFF_EFFECT
baseTable$aa_change <- some93$SNPEFF_AMINO_ACID_CHANGE
UgCl362 <- c(some93$UgCl362 == some93$REF)
UgCl362 <- ifelse(UgCl362== TRUE, "ref", "alt")
baseTable <- cbind(baseTable, UgCl362)
UgCl495 <- c(some93$UgCl495 == some93$REF)
ST93.495 <- ifelse(UgCl495== TRUE, "ref", "alt")
baseTable <- cbind(baseTable, ST93.495)
ST93AB <- baseTable

ST93A.sigALT <- subset(ST93AB, ST93A.ALT==20) #60 variants that are all ALT in ST93A and all REF in ST93B
ST93B.sigALT <-subset(ST93AB, ST93B.ALT==16) #37 variants that are ALT in ST93B and REF in ST93A (2 linked var, two bases apart that are ALT in

ST93A.sigALT.var <- subset(some93, CP %in% ST93A.sigALT$CP)
ST93B.sigALT.var <- subset(some93, CP %in% ST93B.sigALT$CP)

#split into lists for plotting
ST93A.sigALT.var$CHROM <- factor(ST93A.sigALT.var$CHROM,1:14)
ST93B.sigALT.var$CHROM <- factor(ST93B.sigALT.var$CHROM,1:14)

ST93A.sigALT.var.ddn <- split(ST93A.sigALT.var, ST93A.sigALT.var$CHROM)
ST93B.sigALT.var.ddn <- split(ST93B.sigALT.var, ST93B.sigALT.var$CHROM)

table(ST93A.sigALT.var$effect)/60
table(ST93B.sigALT.var$effect)/37

#pdf("manuscript/figures/Figure4A-variantLocationsLine.pdf", width=4, height=7)
tiff(filename = "manuscript/figures/Figure4A-variantLocationsLine.tiff", width = 2, height = 5, units = 'in', res = 300, compression = 'lzw', pointsize = 9)
par(mar=c(1,1,1,1), oma=c(3, 4, 1, 1), fig=c(0, 1, 0, 1), mgp=c(1,0.75,0))
snpLocLine(ST93A.sigALT.var.ddn, ST93B.sigALT.var.ddn)
dev.off()

ST93A.sigALT$clade <- "A"
ST93B.sigALT$clade <- "B"
ST93AB.sig <- rbind(ST93A.sigALT, ST93B.sigALT)
ST93AB.sig <- ST93AB.sig[order(ST93AB.sig$CP),]
ST93AB.sig.var <- subset(some93, CP %in% ST93AB.sig$CP)
ST93AB.sig.var <- ST93AB.sig.var[order(ST93AB.sig.var$CP),]
ST93AB.sig.var$clade <- ST93AB.sig$clade

alias <- c()
for(i in 1:length(ST93AB.sig.var$gene)){
  if(is.na(as.character(ST93AB.sig.var$gene[i]))) alias[i] <- "null"
  else alias[i] <- as.character(subset(genes, name==as.character(ST93AB.sig.var$gene[i]))$alias)
}
ST93AB.sig.var$alias <- alias

description <- c()
for(i in 1:length(ST93AB.sig.var$gene)){
  description[i] <- as.character(subset(genes, name==as.character(ST93AB.sig.var$gene[i]))$description)
}
ST93AB.sig.var$description <- description

ST93AB.sig.var <- ST93AB.sig.var[order(ST93AB.sig.var$gene),]
write.csv(ST93AB.sig.var, "data_out/ST93AB/ST93clade-variants.csv", row.names=FALSE)

ST93AB.sig.var$effect <- as.character(ST93AB.sig.var$effect)

ST93AB.sig.var$effect[ST93AB.sig.var$effect %in% c("CODON_CHANGE_PLUS_CODON_DELETION", "CODON_CHANGE_PLUS_CODON_INSERTION", "CODON_DELETION", "CODON_INSERTION", "FRAME_SHIFT")]<- "indel"
ST93AB.sig.var$effect[ST93AB.sig.var$effect %in% c("START_GAINED", "STOP_GAINED", "NON_SYNONYMOUS_CODING", "STOP_LOST")] <- "nonsynonymous"
ST93AB.sig.var$effect[ST93AB.sig.var$effect %in% c("UTR_5_PRIME", "UPSTREAM")] <- "upstream"
ST93AB.sig.var$effect[ST93AB.sig.var$effect %in% c("UTR_3_PRIME", "DOWNSTREAM")] <- "downstream"
ST93AB.sig.var$effect[ST93AB.sig.var$effect %in% c("SYNONYMOUS_CODING", "INTERGENIC", "INTRON", "SPLICE_SITE_ACCEPTOR", "SPLICE_SITE_DONOR")] <- "synonymous"
ST93AB.sig.var$effect <- factor(ST93AB.sig.var$effect, levels=c("indel", "nonsynonymous", "upstream", "downstream", "synonymous"))

ST93AB.sig.var.effect <- subset(ST93AB.sig.var, effect %in% c("nonsynonymous", "indel"))
write.csv(ST93AB.sig.var.effect, "manuscript/tables/TableS3_ST93clade-variants-effect.csv", row.names=FALSE)

##################################################
#barplots for effect type
##################################################
all93.small <- all93
all93.small$effect <- as.character(all93.small$effect)


all93.small$effect[all93.small$effect %in% c("CODON_CHANGE_PLUS_CODON_DELETION", "CODON_CHANGE_PLUS_CODON_INSERTION", "CODON_DELETION", "CODON_INSERTION", "FRAME_SHIFT")] <- "indel"
all93.small$effect[all93.small$effect %in% c("START_GAINED", "STOP_GAINED", "NON_SYNONYMOUS_CODING", "STOP_LOST")] <- "nonsynonymous"
all93.small$effect[all93.small$effect %in% c("UTR_5_PRIME", "UPSTREAM")] <- "upstream"
all93.small$effect[all93.small$effect %in% c("UTR_3_PRIME", "DOWNSTREAM")] <- "downstream"
all93.small$effect[all93.small$effect %in% c("SYNONYMOUS_CODING", "INTERGENIC", "INTRON", "SPLICE_SITE_ACCEPTOR", "SPLICE_SITE_DONOR")] <- "synonymous"
all93.small$effect <- factor(all93.small$effect, levels=c("indel", "nonsynonymous", "upstream", "downstream", "synonymous"))

ST93A.small <- ST93A.sigALT.var
ST93A.small$effect <- as.character(ST93A.small$effect)
ST93A.small$effect[ST93A.small$effect %in% c("CODON_CHANGE_PLUS_CODON_DELETION", "CODON_CHANGE_PLUS_CODON_INSERTION", "CODON_DELETION", "CODON_INSERTION", "FRAME_SHIFT")] <- "indel"
ST93A.small$effect[ST93A.small$effect %in% c("START_GAINED", "STOP_GAINED", "NON_SYNONYMOUS_CODING", "STOP_LOST")] <- "nonsynonymous"
ST93A.small$effect[ST93A.small$effect %in% c("UTR_5_PRIME", "UPSTREAM")] <- "upstream"
ST93A.small$effect[ST93A.small$effect %in% c("UTR_3_PRIME", "DOWNSTREAM")] <- "downstream"
ST93A.small$effect[ST93A.small$effect %in% c("SYNONYMOUS_CODING", "INTERGENIC", "INTRON", "SPLICE_SITE_ACCEPTOR", "SPLICE_SITE_DONOR")] <- "synonymous"
ST93A.small$effect <- factor(ST93A.small$effect, levels=c("indel", "nonsynonymous", "upstream", "downstream", "synonymous"))

ST93B.small <- ST93B.sigALT.var
ST93B.small$effect <- as.character(ST93B.small$effect)
ST93B.small$effect[ST93B.small$effect %in% c("CODON_CHANGE_PLUS_CODON_DELETION", "CODON_CHANGE_PLUS_CODON_INSERTION", "CODON_DELETION", "CODON_INSERTION", "FRAME_SHIFT")] <- "indel"
ST93B.small$effect[ST93B.small$effect %in% c("START_GAINED", "STOP_GAINED", "NON_SYNONYMOUS_CODING", "STOP_LOST")] <- "nonsynonymous"
ST93B.small$effect[ST93B.small$effect %in% c("UTR_5_PRIME", "UPSTREAM")] <- "upstream"
ST93B.small$effect[ST93B.small$effect %in% c("UTR_3_PRIME", "DOWNSTREAM")] <- "downstream"
ST93B.small$effect[ST93B.small$effect %in% c("SYNONYMOUS_CODING", "INTERGENIC", "INTRON", "SPLICE_SITE_ACCEPTOR", "SPLICE_SITE_DONOR")] <- "synonymous"
ST93B.small$effect <- factor(ST93B.small$effect, levels=c("indel", "nonsynonymous", "upstream", "downstream", "synonymous"))

allEffect <- rbind("all93" = table(all93.small$effect)/length(all93.small[,1]), "ST93A" = table(ST93A.small$effect)/length(ST93A.small[,1]), "ST93B" = table(ST93B.small$effect) /length(ST93B.small[,1]))


#pdf("manuscript/figures/Figure4B_ST93AB-effectVar.pdf", width=6, height=4)
tiff(filename = "manuscript/figures/Figure4B_ST93AB-effectVar.tiff", width = 4, height = 2.5, units = 'in', res = 300, compression = 'lzw', pointsize = 9)
par(mar=c(1,1,1,1), oma=c(3, 4, 1, 1), fig=c(0, 1, 0, 1), mgp=c(1,0.5,0))
#par(mar=c(4, 3, 2, 1), oma=c(2, 2, 1, 1))
m <- barplot(allEffect, beside=TRUE, xaxt="n", yaxt="n", col=c("darkgrey", "purple", "orange"), xaxs = "i", ylim=c(0, 0.8), xlim=c(0, 21), border=NA)
box()
axis(2, las=2, pos=0)
axis(1, at =m[2,], labels=c("indel", "nonsynonymous", "upstream", "downstream", "synonymous"), pos=0, outer=TRUE, cex.axis=0.9)
mtext("fraction of variants", side=2, line=3)
legend("topleft", legend=c("allST93", "ST93A", "ST93B"), pch=22, col=c("darkgrey", "purple", "orange"), pt.bg =c("darkgrey", "purple", "orange"), cex=0.8, inset=0.05)
dev.off()

allEffect_tab <- as.table(rbind(table(ST93A.small$effect), table(ST93B.small$effect)))
chi1 <- chisq.test(allEffect_tab)
