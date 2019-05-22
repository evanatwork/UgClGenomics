#help from:
#https://rgriff23.github.io/2017/05/11/primate-phylogeny-ggtree.html

####################################
#load libraries for use
####################################
library("ape")
# library("Biostrings")
library("ggplot2")
library("ggtree") #bioconductor
library("phangorn")
library(phytools)

####################################
#load tree, tip labels
####################################
tree <- read.raxml("data_in/phylogeny/RAxML_bipartitionsBranchLabels.STall5.bip")

treeOrig <- tree
tree@phylo$tip.label[c(1:5, 20, 21, 37:39, 49:55)] <- c("UgCl021", "UgCl047", "UgCl032", "UgCl029", "UgCl045", "UgCl008", "UgCl001", "UgCl212", "UgCl087", "UgCl107", "UgCl065", "UgCl040", "UgCl037", "UgCl018", "UgCl030", "UgCl057", "UgCl093")

tip_labels <- read.csv("data_in/phylogeny/tip_labels.csv")
#put this back in so that can sort on survival to prune the tree
tip_labels[42, "survival"] <- "survived"
tip_labels[38, "survival"] <- "died"


##########################################
#Plot full tree
##########################################
pdf("manuscript/figures/Figure2A-180718UgCl.pdf", width=5, height=7)
ggtree(tree, ladderize=TRUE)  +
  geom_tiplab(offset=0.001, size=3) +  
  ggplot2::xlim(0, 1)+  
  geom_treescale(x = 0.75, y =2, offset=0.6, width=0.1) +
  geom_point2(aes(label=bootstrap, subset=!is.na(as.numeric(bootstrap)) & bootstrap >50))
#geom_text2(aes(label=bootstrap, subset=!is.na(as.numeric(bootstrap)) & bootstrap >50)), nudge_x = -0.03, nudge_y = 0.2)
dev.off()

##########################################
#Plot just ST93 strains
##########################################
notST93 <- as.character(subset(tip_labels, days=="pre-trial")$strain)
justST93.tree <- drop.tip(treeOrig@phylo, notST93)
justST93.tree$tip.label[29] <- "UgCl212"

tip_labels$strain <- as.character(tip_labels$strain)
tip_labels$strain[tip_labels$strain == "UgCl21"] <- "UgCl212"
sub_labels <- subset(tip_labels, strain %in% justST93.tree$tip.label)[c("strain", "surv.code")]
#x[order(match(x,y))]
sub_labels <- sub_labels[order(match(sub_labels$strain,justST93.tree$tip.label)),]
justST93.tree$survival <- as.character(sub_labels$surv.code)

pdf("manuscript/figures/Figure2B-180718UgClABtree-495.pdf", width=5, height=6)
ggtree(justST93.tree, ladderize=FALSE) +
  geom_tiplab(offset=0.001, size=3)+    
  geom_treescale(x = 0, y =-1, offset=0.6, width=0.01)+
  # theme_tree2(legend.position='right') +
  theme(plot.margin = unit(c(1,3,1,1), "cm")) +
  geom_hilight(node=42, fill="purple", alpha=.4, extend=0.008) +
  geom_hilight(node=61, fill="orange", alpha=.4, extend=0.0058)
dev.off()
