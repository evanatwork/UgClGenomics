library(Hmisc)
lenU <- function(x) length(unique(x))

genes <- read.csv("data_in/general/180709genesList.csv", as.is=TRUE)
strains <- read.csv("data_in/general/180727strainInfo.csv")
strains <- subset(strains, ST=="93" & arm != "pre-COAT")
names(strains)[8] <-"growth"

######################################################
#Read in full variants files (i.e., all variants)
######################################################
dfall <- read.csv("data_in/gatk_processed/UgClSeq_snps.txt", as.is=TRUE)
dfall.IND <- read.csv("data_in/gatk_processed/UgClSeq_indels.txt", as.is=TRUE)
allVar <- rbind(dfall, dfall.IND)
allVar <- subset(allVar, numDat >49) #142376
allVar.table <- table(allVar$SNPEFF_GENE_NAME)

dCSF <- read.csv("data_in/phenotypes/Clean_Compiled_CSF_log2duplicates_d0.csv", as.is=TRUE) #"duplicates", t=0
dCSF$IL1b <- as.numeric(dCSF$IL1b)

dWBC <- read.csv("data_in/phenotypes/COAT_WBC.csv")

######################################################
#Read in variants files (i.e., the possible ST93 variants)
######################################################
var93.dCSF <- read.csv("data_out/variants/var_dCSF.csv", as.is=TRUE)
var93.dCSF$CP <- paste(var93.dCSF$CHROM, var93.dCSF$POS, sep=".")
var93.dCSF$gene <- as.character(var93.dCSF$gene)

var93.dWBC <- read.csv("data_out/variants/var_dWBC.csv", as.is=TRUE)
var93.dWBC$CP <- paste(var93.dWBC$CHROM, var93.dWBC$POS, sep=".")
var93.dWBC$gene <- as.character(var93.dWBC$gene)

var93.dINV  <- read.csv("data_out/variants/var_dINV.csv", as.is=TRUE)
var93.dINV$CP <- paste(var93.dINV$CHROM, var93.dINV$POS, sep=".")
var93.dINV$gene <- as.character(var93.dINV$gene)

######################################################
#Look at all significant survival variants
######################################################
allSigCP <- read.csv("manuscript/tables/TableS6_geneP-allSigCP.csv")

geneP.all <-read.csv("data_out/GWAS/180808genePall.csv")

survival <- subset(geneP.all, Survival < 0.05)[,c(1:6, 26, 42:48)] #12
allSigCP_survival <- subset(allSigCP, gene %in% survival$gene)
#all seven genes with significant survival associations were already included

######################################################
#Bootstrap for significant varints
######################################################
allSig_allPhen <- read.csv("manuscript/tables/TableS7_geneP2-allsig-stats.csv", as.is=TRUE)
allSig_allPhen$CP <- paste(allSig_allPhen$chr, allSig_allPhen$pos, sep=".")
allSig_allPhen$pheno[allSig_allPhen$pheno=="Survival"] <- "survival_days"

j <- 1
traits <- c()
CP <- c()
numFalse.P <- c()
numExtreme <- c()
for(i in unique(allSig_allPhen$CP)){
  subAllSig <- subset(allSig_allPhen, CP==i)
  #eachCP <- geneP.all.sig2[i ,c(1:5,15:40)]
  #sig <- geneP.all.sig2[i, c(14+which(geneP.all.sig2[i,15:40]<0.05))]
  #sub <- subset(someST93.4.eff.imm, CP == eachCP$CP)
  #traits <- append(traits, names(sig))
  traits <- append(traits, subAllSig$pheno)
  CP <- append(CP, rep(subAllSig$CP, nrow(subAllSig)))
  for(l in subAllSig$pheno){
    sub <- subset(subAllSig, pheno == l)
    estims <- c()
    numFalse.P <- c()
    k <- 0
    pvals <- c()
    while(k < 100){
      k <- k+1
      #resp <- sub[15:49]
      if(l %in% names(dCSF)){
        pred <- dCSF[l]
        resp <- subset(var.dCSF, CP == sub$CP)[,1:35] #variant information
      }
      if(l %in% names(dWBC)){
        pred <- dWBC[l]
        resp <- subset(var.dWBC, CP == sub$CP)[,1:37] #variant information
      }
      if(l %in% names(strains)){
        pred <- strains[2:38,l] #phenotypic data
        resp <- subset(var93.dINV, CP == sub$CP)[,1:37] #variant information
      }
      dat <- data.frame(pred, resp = t(resp))
      names(dat) <- c("pred", "resp")
      dat <- subset(dat, !is.na(pred))
      dat <- subset(dat, resp != ".")
      dat$resp <- as.character(dat$resp)
      dat$resp[dat$resp== as.character(subset(var.dCSF, CP == sub$CP)$REF)] <- 1		#REF
      dat$resp[dat$resp== as.character(subset(var.dCSF, CP == sub$CP)$ALT)] <- 0		#ALT
      dat$pred <- as.numeric(dat$pred)
      dat$resp <- as.numeric(dat$resp)
      dat.real <- dat
      dat$resp <- sample(dat$resp)
      t2 <- glm(resp~pred, data=dat, family=binomial)
      estims[k] <- summary(t2)$coefficients[2,1]
      pvals[k] <- summary(t2)$coefficients[2,4]
    }
    #numFalse.P <- append(numFalse.P, table(pvals > 0.05)["FALSE"])
    #print(table(pvals < 0.05)["FALSE"])
    t2.real <- glm(resp~pred, data=dat.real, family=binomial)
    estim.real <- summary(t2.real)$coefficients[2,1]
    numHigh <- length(subset(estims, estims > estim.real))
    numLow <- length(subset(estims, estims < estim.real))
    numExtreme <- append(numExtreme, min(numHigh, numLow))
    #plot(estims)
    #abline(h=summary(t2.real)$coefficients[2,1], col="red")
    j <- j+1
  }
}

hist(numExtreme, breaks = 40, xlab = "randomized model is more extreme then measured data", main="", yaxt="n", xaxt="n")
axis(2, las=2, pos=-0.5)
axis(1, pos=0)
abline(v=5, col="red")s

allSig_allPhen$numExtreme <- numExtreme
write.csv(allSig_allPhen, "data_out/GWAS/GWAS_bootstrap.csv")

allSig_allPhen_remove <- subset(allSig_allPhen, numExtreme > 5)

plot(allSig_allPhen$nonSig, allSig_allPhen$numExtreme)

j <- 1
traits <- c()
CP <- c()
numFalse.P <- c()
estim.smaller <- c()
estim.larger <- c()
nonSig <- c()
for(i in unique(allSig_allPhen$CP)){
  subAllSig <- subset(allSig_allPhen, CP==i)
  traits <- append(traits, subAllSig$pheno)
  CP <- append(CP, rep(subAllSig$CP, nrow(subAllSig)))
  for(l in subAllSig$pheno){
    sub <- subset(subAllSig, pheno == l)
    estims <- c()
    numFalse.P <- c()
    k <- 0
    pvals <- c()
    if(l %in% names(dCSF)){
        pred <- dCSF[l]
        resp <- subset(var.dCSF, CP == sub$CP)[,1:35] #variant information
      }
      if(l %in% names(dWBC)){
        pred <- dWBC[l]
        resp <- subset(var.dWBC, CP == sub$CP)[,1:37] #variant information
      }
      if(l %in% names(strains)){
        pred <- strains[2:38,l] #phenotypic data
        resp <- subset(var93.dINV, CP == sub$CP)[,1:37] #variant information
      }
      dat <- data.frame(pred, resp = t(resp))
      names(dat) <- c("pred", "resp")
      dat <- subset(dat, !is.na(pred))
      dat <- subset(dat, resp != ".")
      dat$resp <- as.character(dat$resp)
      dat$resp[dat$resp== as.character(subset(var.dCSF, CP == sub$CP)$REF)] <- 1		#REF
      dat$resp[dat$resp== as.character(subset(var.dCSF, CP == sub$CP)$ALT)] <- 0		#ALT
      dat$pred <- as.numeric(dat$pred)
      dat$resp <- as.numeric(dat$resp)
      dat.real <- dat
      estims <- c()
      pvals <- c()
      for(k in 1:nrow(dat)){
        temp <- dat[-k,]
        t2 <- glm(resp~pred, data=temp, family=binomial)
        summary(t2)
        #estims[k] <- summary(t2)$coefficients[2,1]
        pvals[k] <- summary(t2)$coefficients[2,4]
      }
      nonSig <- append(nonSig, length(subset(pvals, pvals > 0.05)))
    }
    #numFalse.P <- append(numFalse.P, table(pvals > 0.05)["FALSE"])
    #print(table(pvals < 0.05)["FALSE"])
    # t2.real <- glm(resp~pred, data=dat.real, family=binomial)
    # estim.smaller <- append(estim.smaller, table(summary(t2.real)$coefficients[2,1] < estims)["FALSE"])
    # estim.larger <- append(estim.larger, table(summary(t2.real)$coefficients[2,1] > estims)["TRUE"])
    # j <- j+1
  }
}

allSig_allPhen$nonSig <- nonSig

j <- 1
traits <- c()
CP <- c()
deltas <- c()
for(i in unique(allSig_allPhen$CP)){
  subAllSig <- subset(allSig_allPhen, CP==i)
  traits <- append(traits, subAllSig$pheno)
  CP <- append(CP, rep(subAllSig$CP, nrow(subAllSig)))
  for(l in subAllSig$pheno){
    sub <- subset(subAllSig, pheno == l)
    estims <- c()
    numFalse.P <- c()
    k <- 0
    pvals <- c()
    if(l %in% names(dCSF)){
      pred <- dCSF[l]
      resp <- subset(var.dCSF, CP == sub$CP)[,1:35] #variant information
    }
    if(l %in% names(dWBC)){
      pred <- dWBC[l]
      resp <- subset(var.dWBC, CP == sub$CP)[,1:37] #variant information
    }
    if(l %in% names(strains)){
      pred <- strains[2:38,l] #phenotypic data
      resp <- subset(var93.dINV, CP == sub$CP)[,1:37] #variant information
    }
    dat <- data.frame(pred, resp = t(resp))
    names(dat) <- c("pred", "resp")
    dat <- subset(dat, !is.na(pred))
    dat <- subset(dat, resp != ".")
    dat$resp <- as.character(dat$resp)
    dat$resp[dat$resp== as.character(subset(var.dCSF, CP == sub$CP)$REF)] <- 1		#REF
    dat$resp[dat$resp== as.character(subset(var.dCSF, CP == sub$CP)$ALT)] <- 0		#ALT
    dat$pred <- as.numeric(dat$pred)
    dat$resp <- as.numeric(dat$resp)

    fit <- glm(resp~pred, data=dat, family=binomial)
    out <- NULL
    for (m in 1:nrow(dat))
      out[m] <- predict(update(fit, data = dat[-m,]), newdata = dat[m,], type = "response")

    deltas <- append(deltas, boot::cv.glm(dat, fit)$delta[1])
    ## 0.1497903 = mtcars
    # 0.28
    #mean((dat$resp - round(out))^2)
    ## 0.1875 - mtcars
    # 0.515
    #mean((dat$resp - out)^2)
    ## 0.1497903- mtcars
    # 0.28
    }
}
allSig_allPhen$deltas <- deltas

write.csv(allSig_allPhen, "data_out/GWAS/GWAS_dropOne.csv")

plot(jitter(as.numeric(as.factor(allSig_allPhen$pheno))), allSig_allPhen$nonSig)

plot(as.factor(allSig_allPhen$pheno), allSig_allPhen$nonSig)

plot(allSig_allPhen$nonSig, allSig_allPhen$estim.larger, ylab="bootstrap estimates more extreme", xlab= "nonsig with drop one strain")
abline(h=95, lty=2, col="red")
abline(h=5, lty=2, col="red")
abline(v=24.5, lty=2, col="red")

plot(allSig_allPhen$nonSig, allSig_allPhen$deltas)
abline(v=24.5, lty=2, col="red")



#############################
#https://stats.stackexchange.com/questions/250809/confusion-about-cv-glm-in-r
fit <- glm(vs ~ mpg, data = mtcars, family = binomial)
out <- NULL
for (i in 1:nrow(mtcars))
  out[i] <- predict(update(fit, data = mtcars[-i,]), newdata = mtcars[i,], type = "response")

boot::cv.glm(mtcars, fit)$delta[1]
## 0.1497903
mean((mtcars$vs - round(out))^2)
## 0.1875
mean((mtcars$vs - out)^2)
## 0.1497903
