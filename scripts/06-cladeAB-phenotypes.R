library(here)
library(Hmisc)

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
#Read in genotype files
###########################
ST93AB.sig <-read.csv("tables_intermediate/ST93AB/ST93clade-variants.csv")

###########################
#Read in phenotype files
###########################
dCSF.dup.0 <- read.csv("data_in/phenotypes/Clean_Compiled_CSF_log2duplicates_d0.csv", as.is=TRUE)
dCSF.dup.0 <- subset(dCSF.dup.0, line!="468")
dCSF.dup.0 <- subset(dCSF.dup.0, line!="362")
dCSF.dup.0 <- subset(dCSF.dup.0, line!="495")


dWBC <- read.csv("data_in/phenotypes/COAT_WBC.csv")
dWBC <- subset(dWBC, line!="362")
dWBC <- subset(dWBC, line!="495")

###########################
#T-test for CSF results
############################
ttest.CSF.dup<- data.frame("Phenotype" = c("temp"), "W" = 1, "pvalue" = 1, "A-B" =1)
ttest.CSF.dup$Phenotype <- as.character(ttest.CSF.dup$Phenotype)

#Not enough data MIP1a
k<-0
CSFcols <- c(5:22, 24, 25)
for(i in CSFcols){
  k <- k+1
  test <- wilcox.test(subset(dCSF.dup.0, ST93=="A")[,i], subset(dCSF.dup.0, ST93=="B")[,i])
  test2 <- t.test(subset(dCSF.dup.0, ST93=="A")[,i], subset(dCSF.dup.0, ST93=="B")[,i])
  ttest.CSF.dup[k,] <- c(as.character(names(dCSF.dup.0)[i]), round(test$statistic, 2),  round(test$p.value, 3), round(c(test2$estimate[1]-test2$estimate[2]), 2))
}
names(ttest.CSF.dup) <- c("phenotype", "W", "pvalue", "A-B")
#write.csv(ttest.CSF.dup, "manuscript/tables/TableS3_180727CSFdup-clade_ttest.csv", row.names=FALSE)

ttest.WBC<- data.frame("Phenotype" = c("temp"), "W" = 1, "pvalue" = 1, "A-B" =1)
ttest.WBC$Phenotype <- as.character(ttest.WBC$Phenotype)

WBCcols <- c(5:9)
k<-0
for(i in WBCcols){
  k <- k+1
  test <- wilcox.test(subset(dWBC, ST93=="A")[,i], subset(dWBC, ST93=="B")[,i])
  test2 <- t.test(subset(dWBC, ST93=="A")[,i], subset(dWBC, ST93=="B")[,i])
  ttest.WBC[k,] <- c(names(dWBC)[i], round(test$statistic, 2),  round(test$p.value, 3), round(c(test2$estimate[1]-test2$estimate[2]), 2))
}
names(ttest.WBC) <- c("phenotype", "W", "pvalue", "A-B")
#write.csv(ttest.WBC, "tables/ST93AB/180727WBC-clade_ttest.csv", row.names=FALSE)

ttest.invitro<- data.frame("Phenotype" = c("temp"), "W" = 1, "pvalue" = 1, "A-B" =1)
ttest.invitro$Phenotype <- as.character(ttest.invitro$Phenotype)

invitrocols <- c(5:11)
k<-0
for(i in invitrocols){
  k <- k+1
  test <- wilcox.test(subset(strains, ST93=="A")[,i], subset(strains, ST93=="B")[,i])
  test2 <- t.test(subset(strains, ST93=="A")[,i], subset(strains, ST93=="B")[,i])
  ttest.invitro[k,] <- c(names(strains)[i], round(test$statistic, 2), round(test$p.value, 3), round(c(test2$estimate[1]-test2$estimate[2]), 2))
}
names(ttest.invitro) <- c("phenotype", "W", "pvalue", "A-B")
#write.csv(ttest.invitro, "tables/ST93AB/180727invitro-clade_ttest.csv", row.names=FALSE)

stats.ttest <- rbind(ttest.CSF.dup, ttest.WBC, ttest.invitro)
write.csv(stats.ttest, "manuscript/tables/TableS4_ST93ABallPheno-ttest.csv", row.names=FALSE)
#################################
#Figures
#################################
pdf("manuscript/figures/Figure4C_sigPhenoAB.pdf", width=4, height=4.5)
par(mfrow=c(2, 1), mar=c(1, 3, 1, 1), oma=c(3, 1, 1, 1))
plot(jitter(as.numeric(as.factor(subset(strains, ST93 %in% c("A", "B"))$ST93))), as.numeric(subset(strains, ST93 %in% c("A", "B"))$uptake), xaxt="n", yaxt="n", xlim=c(0.5, 2.5), ylim=c(0, 1.6), main="", cex=1.2)
arrows(0.8, median(subset(strains, ST93== "A")$uptake, na.rm=TRUE), 1.2, median(subset(strains, ST93== "A")$uptake, na.rm=TRUE), lwd=3,
col="purple", length=0)
arrows(1.8,median(subset(strains, ST93== "B")$uptake, na.rm=TRUE), 2.2, median(subset(strains, ST93== "B")$uptake, na.rm=TRUE), lwd=3, col="goldenrod", length=0)
axis(2, las=2)
axis(1, at=1:2, labels=FALSE)
mtext("uptake", side=3, adj=0.01, font=2)

plot(jitter(as.numeric(as.factor(subset(dCSF.dup.0, ST93 %in% c("A", "B"))$ST93))), as.numeric(subset(dCSF.dup.0,  ST93 %in% c("A", "B"))$IL2), xaxt="n", yaxt="n", ylim=c(-5, 6), main="", xlim=c(0.5, 2.5), ylab="", xlab="")
axis(2, las=2)
arrows(0.8, median(subset(dCSF.dup.0, ST93== "A")$IL2, na.rm=TRUE), 1.2, median(subset(dCSF.dup.0, ST93== "A")$IL2, na.rm=TRUE), lwd=3,
col="purple", length=0)
arrows(1.8, median(subset(dCSF.dup.0, ST93== "B")$IL2, na.rm=TRUE), 2.2, median(subset(dCSF.dup.0, ST93== "B")$IL2, na.rm=TRUE), lwd=3, col="goldenrod", length=0)
axis(1, at=1:2, labels=c("ST93A", "ST93B"))
mtext("IL2", side=3, adj=0.01, font=2)
dev.off()

yrange <- data.frame(min = c(-5, -5, -5, -3, 0, -6, 5, 0, -2, -2, -1, 2, 7.5, -2, 5, 3, 0, 2, 0, 0, 0, 0, 0, -0.38, 0, 0, 0, 0, 0, 0, 0), max = c(4, 6, 5, 5, 14, 5, 13, 7, 7, 10, 6, 9, 10, 10, 15, 10, 7, 9, 35000, 600, 1000000, 120, 400, -0.22, 1.5, 1.6, 7500, 0.4, 70, 1.2, 13))
labels <-c(names(dCSF.dup.0)[CSFcols], "HIV-RNA", "CD4", "CSF WBC", "EFA", "uptake", "adherence", "chitin", "absolute growth", "fluconazole MIC", "amphotericin B MIC", "sertraline MIC")
pdf("manuscript/figures/FigureS1_allPheno-cladeComparison.pdf", width=9, height=10)
  k <- 0
  par(mfrow=c(8, 4), mar=c(1, 3, 1, 1), oma=c(3, 2, 1, 1))
  for(i in CSFcols){
    k <- k+1
    plot(jitter(as.numeric(as.factor(subset(dCSF.dup.0, ST93 %in% c("A", "B"))$ST93))), as.numeric(subset(dCSF.dup.0,  ST93 %in% c("A", "B"))[,i]), xaxt="n", yaxt="n", ylim=c(yrange[k,1], yrange[k,2]), xlim=c(0.5, 2.5), ylab="")
    axis(2, las=2)
    arrows(0.8, median(subset(dCSF.dup.0, ST93== "A")[,i], na.rm=TRUE), 1.2, median(subset(dCSF.dup.0, ST93== "A")[,i], na.rm=TRUE), lwd=3,
    col="purple", length=0)
    arrows(1.8, median(subset(dCSF.dup.0, ST93== "B")[,i], na.rm=TRUE), 2.2, median(subset(dCSF.dup.0, ST93== "B")[,i], na.rm=TRUE), lwd=3, col="goldenrod", length=0)
    axis(1, labels=FALSE, at=c(1, 2))
    mtext(paste0(labels[k], ", p = ", stats.ttest[k, "pvalue"]), adj=0.01, cex=0.9)
  }
  for(j in c(6:9)){
    k <- k+1
    plot(jitter(as.numeric(as.factor(dWBC$ST93))), as.numeric(dWBC[,j]), xaxt="n", yaxt="n", ylim=c(yrange[k,1], yrange[k,2]), xlim=c(0.5, 2.5))
    axis(2, las=2)
    arrows(0.8, median(subset(dWBC, ST93== "A")[,j], na.rm=TRUE), 1.2, median(subset(dWBC, ST93== "A")[,j], na.rm=TRUE), lwd=3,
    col="purple", length=0)
    arrows(1.8, median(subset(dWBC, ST93== "B")[,j], na.rm=TRUE), 2.2, median(subset(dWBC, ST93== "B")[,j], na.rm=TRUE), lwd=3, col="goldenrod", length=0)
    axis(1, labels=FALSE, at=c(1, 2))
    mtext(paste0(labels[k], ", p = ", stats.ttest[k, "pvalue"]), adj=0.01, cex=0.9)
}
for(l in invitrocols){
  k <- k+1
  plot(jitter(as.numeric(as.factor(subset(strains, ST93 %in% c("A", "B"))$ST93))), as.numeric(subset(strains, ST93 %in% c("A", "B"))[,l]), xaxt="n", yaxt="n", xlim=c(0.5, 2.5), ylim=c(yrange[k,1], yrange[k,2]),  cex=1.2)
  arrows(0.8, median(subset(strains, ST93== "A")[,l], na.rm=TRUE), 1.2, median(subset(strains, ST93== "A")[,l], na.rm=TRUE), lwd=3,
  col="purple", length=0)
  arrows(1.8, median(subset(strains, ST93== "B")[,l], na.rm=TRUE), 2.2, median(subset(strains, ST93== "B")[,l], na.rm=TRUE), lwd=3, col="goldenrod", length=0)
  mtext(paste0(labels[k], ", p = ", stats.ttest[k, "pvalue"]), adj=0.01, cex=0.9)
  axis(2, las=2)
  if(k > 27) axis(1, at= c(1, 2), labels=c("clade-A", "clade-B"), las=1, cex.axis=1.2)
  else axis(1, labels=FALSE, at=c(1, 2))
 }
  dev.off()
