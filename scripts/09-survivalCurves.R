library(survival)
library(dplyr)
#library(OIsurv) # Aumatically loads KMsurv
#library(ranger)
library(ggplot2)
library(survminer)
library(broom)

#https://www.datacamp.com/community/tutorials/survival-analysis-R

#all data for plotting
surv <- read.csv("data_in/KO/180910deletionSurvival.csv")

#data for analysis
survD <- read.csv("data_in/KO/180907KO_daySurvival.csv", as.is=TRUE)

#split by experiment because have to compare to KN99 within experiment
###EXPERIMENT 1
survD.1 <- subset(survD, exper=="1")

#create a set of lists that each contain the control strain KN99 and one other strain
ddn.1 <- split(survD.1, survD.1$strain)

pltList <- list()
fit <- list()
sdiff1 <- list()
for (i in 1:14){
  df <- rbind(ddn.1[[i]], ddn.1[[15]]) #KN99 is in 15
  surv_object <- Surv(time = df$daySurv, event = df$censored)
  fit <- survfit(surv_object~df$strain)
  sdiff1[[i]] <- survdiff(surv_object~df$strain)
  pltList[[i]] <- ggsurvplot(fit, data=df, pval=TRUE, pval.method=TRUE)
}



#this sort of works
 arrange_ggsurvplots(pltList, ncol=4, nrow=4)

#just pull the pvalues from the individual plots manuyally and then make my own figures
pval1 <- c(0.0027, 0.82, 0.08, 0.095, 0.0016, 0.079, 0.43, 0.79, 0.72, 0.044, 0.0027, 0.016, 0.82, 0.18)

strain1 <- c("2176", "363","4373", "4535", "4922", "5662", "5663", "5913", "6169", "6332", "6574", "6704", "6876", "7837")

###EXPERIMENT 2
survD.2 <- subset(survD, exper=="2")

#create a set of lists that each contain the control strain KN99 and one other strain
ddn.2 <- split(survD.2, survD.2$strain)

pltList2 <- list()
sdiff2 <- list()
for (i in 1:4){
  df <- rbind(ddn.2[[i]], ddn.2[[5]]) #KN99 is in 45
  surv_object <- Surv(time = df$daySurv, event = df$censored)
  fit <- survfit(surv_object~df$strain)
  sdiff2[[i]] <- survdiff(surv_object~df$strain)
  pltList2[[i]] <- ggsurvplot(fit, data=df, pval=TRUE)
}


pltList2[[1]]

#just pull the pvalues from the individual plots manuyally and then make my own figures
pval2 <- c(0.77, 0.31, 0.0082, 0.31)
strain2 <- c("5937", "6490", "6986","ITR4comp")

###EXPERIMENT =3
survD.3 <- subset(survD, exper=="3")

surv_object3 <- Surv(time = survD.3$daySurv, event = survD.3$censored)
fit3 <- survfit(surv_object~survD.3$strain)
p3 <- ggsurvplot(fit3, data=survD.3, pval=TRUE)
sdiff3 <- survdiff(surv_object3~survD.3$strain)

#just pull the pvalues from the individual plots manuyally and then make my own figures
pval3 <- c(0.31)
strain3 <- c("7703")

###EXPERIMENT 4
survD.4 <- subset(survD, exper=="4")

#create a set of lists that each contain the control strain KN99 and one other strain
ddn.4 <- split(survD.4, survD.4$strain)

pltList4 <- list()
sdiff4 <- list()
for (i in 1:2){
  df <- rbind(ddn.4[[i]], ddn.4[[3]]) #KN99 is in 3
  surv_object <- Surv(time = df$daySurv, event = df$censored)
  fit <- survfit(surv_object~df$strain)
  pltList4[[i]] <- ggsurvplot(fit, data=df, pval=TRUE)
  sdiff4[[i]] <- survdiff(surv_object~df$strain)
}

pltList4[[1]]

#just pull the pvalues from the individual plots manuyally and then make my own figures
pval4 <- c(0.47, 0.013)
strain4 <- c("ITR4comp", "5662")

#bring all experiments together
stats <- rbind(glance(sdiff1[[1]]), glance(sdiff1[[2]]), glance(sdiff1[[3]]), glance(sdiff1[[4]]), glance(sdiff1[[5]]), glance(sdiff1[[6]]), glance(sdiff1[[7]]), glance(sdiff1[[8]]), glance(sdiff1[[9]]), glance(sdiff1[[10]]), glance(sdiff1[[11]]), glance(sdiff1[[12]]), glance(sdiff1[[13]]), glance(sdiff1[[14]]), glance(sdiff2[[1]]), glance(sdiff2[[2]]), glance(sdiff2[[3]]), glance(sdiff2[[4]]), glance(sdiff3), glance(sdiff4[[1]]), glance(sdiff4[[2]]))

df.ps <- data.frame(strain=c(strain1, strain2, strain3, strain4), pval = c(pval1, pval2, pval3, pval4), exper = c(rep(1, length(strain1)), rep(2, length(strain2)), rep(3, length(strain3)), rep(4, length(strain4))), cbind(stats))
df.sig <- subset(df.ps, pval < 0.05)
df.nsig  <- subset(df.ps, pval > 0.05)

#Table 4
write.csv(df.ps, "manuscript/tables/Table4_1902KOsurvival-stats.csv", row.names=FALSE)
############
#Figures
############

surv.sub <- subset(surv, strain %in% df.sig$strain)

addPoints <- function(i, colour, j = 1, lty=1, wid=1.5){
  points(subset(surv, strain==i & exper == j)$days, subset(surv, strain==i  & exper == j)$numSurv, col=colour, type="l", lwd=wid, lty=lty)
  }

#Sig
pdf("manuscript/figures/Figure6_180910KOsurvival.pdf", width=8, height=4)
par(mar=c(1, 3, 1, 1), oma=c(3, 3, 1, 1), mgp=c(1.5, 0.6, 0))
plot(subset(surv, strain=="KN99" & exper==1 )$days, subset(surv, strain=="KN99" & exper==1 )$numSurv, col="black", type="l", lwd=2, lty=3, yaxt="n", ylab="", xlab="", xaxt="n")
axis(2, las=2, cex.axis=1)
axis(1, cex.axis=1)
addPoints("2176", "red4")
addPoints("6574", "red")
addPoints("6332", "goldenrod")
addPoints("4922", "navyblue")
addPoints("6986", "cornflowerblue", 2)
addPoints("KN99", grey(0.2), 2, lty=2, wid=2)
legend("bottomleft", legend=c(expression(paste("KN99", alpha, " (E1)")), expression(paste("KN99", alpha, " (E2)")), "CNAG_2176 (E1)", "CNAG_6574 (E1)", "CNAG_6332 (E1)", "CNAG_6968 (E2)", "CNAG_4922 (E1)"), col=c("black", grey(0.2), "red4", "red", "goldenrod",  "cornflowerblue", "navyblue"), lwd=c(2, 2, 1, 1, 1, 1, 1), lty=c(3, 2, 1, 1, 1, 1, 1), cex=0.8, bty="n")
mtext("number surviving mice", side=2, line=-1, cex=1.2, outer=TRUE)
mtext("days post infection", side=1, line=2, cex=1.2)
dev.off()

itr4 <- subset(surv, exper=="ITR")
#ITR4
pdf("manuscript/figures/Figure7A_180910ITR4survival.pdf", width=5, height=4)
plot(subset(itr4, strain=="KN99")$days, subset(itr4, strain=="KN99")$numSurv, col="black", type="l", lwd=2, lty=3, yaxt="n", ylab="", xlab="", xaxt="n", xlim=c(0, 45))
points(subset(itr4, strain=="ITR4comp")$days, subset(itr4, strain=="ITR4comp")$numSurv, col=grey(0.4), type="l", lwd=2.5)
points(subset(itr4, strain=="itr4D")$days, subset(itr4, strain=="itr4D")$numSurv, col="dodgerblue", type="l", lwd=2.5)
legend("bottomleft", legend=c(expression(paste("KN99", alpha)), expression(paste("CNAG_5662 (",italic(itr4), Delta, ")")), expression(paste(italic(itr4), Delta, ":", italic(ITR4)))), col=c("black",  "dodgerblue", grey(0.4)), lwd=c(2, 2.5, 2.5), cex=0.8, lty=c(3, 1, 1), bty="n")
axis(2, las=2, cex.axis=1)
axis(1, cex.axis=1)
mtext("Number of Mice", side=2, line=-1, cex=1.2, outer=TRUE)
mtext("Days post Infection", side=1, line=2, cex=1.2)
dev.off()


#non-sig
surv$strain <- as.character(surv$strain)
labels <- c("CNAG_00363", "CNAG_04373", "CNAG_04535", "CNAG_05663", "CNAG_05913", "CNAG_06169", "CNAG_06876", "CNAG_07837", "CNAG_05937", "CNAG_06490", "CNAG_07703")
pdf("manusript/figures/FigureS3_180912KOsurvival-nonSig.pdf", width=8, height=6)
par(mfrow=c(3, 4), mar=c(1,1, 1, 1), oma=c(3, 3, 1, 1), mgp=c(1.5, 0.6, 0))
k <- 0
for(i in c("363", "4373", "4535", "5663", "5913", "6169", "6876", "7837")){
  k <- k+1
  temp <- subset(surv, strain==i)
  plot(temp$days, temp$numSurv, col=grey(0.4), type="l", lwd=2, lty=1, yaxt="n", ylab="", xlab="", xaxt="n", xlim=c(0, 44))
  points(subset(surv, strain=="KN99" & exper==1)$days, subset(surv, strain=="KN99" & exper==1)$numSurv, col="black", type="l", lwd=2, lty=3)
  mtext(paste0(labels[k], " (E1)"), side=3, adj=0.01, cex=0.8, line=0.25)
  if (k %% 4 == 1) axis(2, las=2, cex.axis=1)
  else (axis(2, labels=FALSE))
  if(k > 7) axis(1)
  else (axis(1, labels=FALSE))
}
for(i in c("5937", "6490")){
  k <- k+1
  temp <- subset(surv, strain==i)
  plot(temp$days, temp$numSurv, col=grey(0.4), type="l", lwd=2, lty=1, yaxt="n", ylab="", xlab="", xaxt="n", xlim=c(0, 44))
  points(subset(surv, strain=="KN99" & exper==2)$days, subset(surv, strain=="KN99" & exper==2)$numSurv, col="black", type="l", lwd=2, lty=3)
  mtext(paste0(labels[k], " (E2)"), side=3, adj=0.01, cex=0.8, line=0.25)
  if (k %% 4 == 1) axis(2, las=2, cex.axis=1)
  else (axis(2, labels=FALSE))
  if(k > 7) axis(1)
  else (axis(1, labels=FALSE))
}
for(i in c("7703")){
  k <- k+1
  temp <- subset(surv, strain==i)
  plot(temp$days, temp$numSurv, col=grey(0.4), type="l", lwd=2, lty=1, yaxt="n", ylab="", xlab="", xaxt="n", xlim=c(0, 44))
  points(subset(surv, strain=="KN99" & exper==3)$days, subset(surv, strain=="KN99" & exper==3)$numSurv, col="black", type="l", lwd=2, lty=3)
  mtext(paste0(labels[k], " (E3)"), side=3, adj=0.01, cex=0.8, line=0.25)
  if (k %% 4 == 1) axis(2, las=2, cex.axis=1)
  else (axis(2, labels=FALSE))
  if(k > 7) axis(1)
  else (axis(1, labels=FALSE))
}
mtext("number surviving mice", side=2, line=1, cex=1.2, outer=TRUE)
mtext("days post infection", side=1, line=1.75, cex=1.2, outer=TRUE)
plot(1:1, ann=FALSE, bty="n", xaxt="n", yaxt="n", type="n")
legend("center", legend=c(expression(paste("KN99", alpha)), "knock-out strain"), lty=c(2, 1))
dev.off()

#find aliases
