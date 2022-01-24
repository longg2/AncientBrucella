library(dplyr)
library(ggplot2)
library(ggtree)
library(phangorn)
library(treeio)
library(ggpubr)
library(ggstance)
library(ggnewscale)
library(phytools)
library(ape)
library(gtools)

# For the colours
colour <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6',
'#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3',
'#808000', '#ffd8b1', '#000075', '#808080', '#f0f0ff', '#000000')

ann_colors <- list(ST = c("5" = colour[10], "7" = colour[15], "8" = colour[3], "9" = colour[9], "10" = colour[6], "11" = colour[17], "12" = colour[8], "39" = colour[13], "40" = colour[2], "41" = colour[16], "42" = colour[7], "43" = colour[4], "71" = colour[11], "88" = colour[12], "102" = colour[19], "Ancient" = colour[1], "NF" = colour[20], "Reference" = colour[22], "Kay" = colour[5]))
##############

# Basic Tree
tree <- read.newick(file = "ST11KayIncluded.nwk")
tree <- phytools::reroot(tree, interactive = T)
Tiplabels <- tree$tip.label

##################################################
# Let's colour these based on their PhyloGroups ##
##################################################
STData <- read.delim("MLSTResults.txt")[,1:2]
rownames(STData) <- STData$Sample
STData["JessSample",] <- list(Sample = "JessSample", ST = "Ancient")
STData["Reference",] <- list(Sample = "Reference", ST = "Reference")
STData["KayBMel",] <- list(Sample = "KayBMel", ST = "Kay")

STData <- STData[tree$tip.label,]
STData <- STData %>% mutate(ST = gsub("[^[:alnum:] ]","", ST))

# Let's filter the ann_colors list so that only the STs present in the phylogeny are used
ann_colors$ST <- ann_colors$ST[unique(STData$ST)]

# Going to try replacing the labels with the Phylogroups
p1 <- ggtree(tree, right = T) %<+% STData +
       	geom_nodepoint(size = 1.5, colour = ifelse(as.numeric(tree$node.label) >= 90, "black", ifelse(as.numeric(tree$node.label) >= 50, "grey",NA)),
		     shape = "square") +
	geom_tippoint(mapping = aes(colour = ST), size = 1.5) +
	#geom_tiplab(align = T) + 
	geom_treescale(linesize = 1, offset = 2) + theme_tree() + geom_rootedge(rootedge = 0.00005) +
	scale_color_manual(values = ann_colors$ST, name = "Sequence Type") +
	theme(legend.position = "bottom")# + guides(colour = guide_legend(nrow = 1))
p1
ggsave(p1, file = "BMelCoreST11SNPPhyloKay.pdf", width = 8, height = 6)

#ggsave(p2, file = "~/Documents/University/EcoliPaperV2/AdditionalFiles/PartsofFigures/PathPhylo.pdf", width = 12, height = 9)
load("FullPhylo.RData")
load("ST4995Only.RData")
p23 <- ggarrange(p1,p3, ncol = 2, labels = c("B", "C"), align = "hv")

p123 <- ggarrange(p4, p23, nrow =2, labels = c("A",""), align = "hv", heights = c(4,3), common.legend = F, legend = "bottom")
ggsave(plot = p123, "~/Documents/University/EcoliPaperV2/Figures/Figure3ML.pdf", width = 12, height = 8)

#write.tree(tree, file = "MidpointTreemerTree.nwk")

########################
### The Tempest Plot ###
########################

# Now to get the Tempest Plot
dat <- as_tibble(read.delim("PhylogenyNoContigsTempest.tab")) %>% filter(date > 0)
# Getting the Phylogroups plotted

STData <- STData[match(dat$tip, STData$Sample),] # Getting them in the Right Order
dat$ST <- STData$ST

tempest <- dat %>% ggplot(aes(x = date, y  = distance)) + geom_smooth(method = "lm", show.legend = F, colour = "black") +
	geom_point(alpha = 0.75, aes(colour = ST), show.legend = T) + theme_bw() + ylab("Root to Tip Divergence") + xlab("Year") +
	scale_color_manual(values = ann_colors$ST, name = "Sequence Type")

model <- lm(distance ~date, data = dat)
summary(model)
r2 <- round(summary(model)$adj.r.squared,3)
tmp <- summary(model)$fstatistic
pval <- round(pf(tmp[1],tmp[2],tmp[3], lower.tail = F),3)

tempest <- tempest + annotate(geom = "text", x = 1700, y = 0.0018,label = bquote(R[adj]^2 == .(r2))) +
	annotate(geom = "text", x = 1700, y = 0.0017, label = bquote(P == .(pval))) +
	scale_y_continuous(breaks = scales:::pretty_breaks(n = 10)) + theme(legend.position = "bottom") + coord_cartesian(ylim = c(0,0.002))
tempest
ggsave("FullPhyloTempestNOCONTIGS.pdf", width = 6, height = 4)

load("TempestAll.RData")
load("TempestST4995.RData")

ggarrange(tempestAll, tempest,tempestST4995ONLY, ncol = 1, align = "hv", labels = "AUTO")
ggsave("~/Documents/University/EcoliPaperV2/Figures/TempestPlots.pdf", width = 9, height = 12)

# Now to make the distance matrix
d <- cophenetic.phylo(tree)

subsetKeep <- c("AncientEcoli","ATCC11229", "STEC388","035-003-ccl", "IAI1")
d[subsetKeep,subsetKeep] %>% xtable:::xtable() %>% print(file = "~/Documents/University/EcoliPaperV2/AdditionalFiles/DistMatrix.tex")


