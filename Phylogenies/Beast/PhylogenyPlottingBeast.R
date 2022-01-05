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

ann_colors <- list(ST = c("5" = colour[10], "7" = colour[15], "8" = colour[3], "9" = colour[9], "10" = colour[6], "11" = colour[17], "12" = colour[8], "39" = colour[13], "40" = colour[2], "41" = colour[16], "42" = colour[7], "43" = colour[4], "71" = colour[5], "88" = colour[12], "102" = colour[19], "Ancient" = colour[1], "NF" = colour[20], "Reference" = colour[22]))

# Basic Tree
tree <- read.beast(file = "BeastResults/FinalST11BeastTree.nexus")
Tiplabels <- tree@phylo$tip.label

##################################################
# Let's colour these based on their PhyloGroups ##
##################################################
STData <- read.delim("../MLSTResults.txt")[,1:2]
rownames(STData) <- STData$Sample
STData["JessSample",] <- list(Sample = "JessSample", ST = "Ancient")
STData["Reference",] <- list(Sample = "Reference", ST = "Reference")

STData<- STData[Tiplabels,]

# Let's filter the ann_colors list so that only the STs present in the phylogeny are used
ann_colors$ST <- ann_colors$ST[unique(STData$ST)]

# Going to try replacing the labels with the Phylogroups
p1 <- ggtree(tree, right = T, mrsd = "2017-01-01") %<+% STData +
	geom_tippoint(mapping = aes(colour = ST), size = 1.5) +
	geom_range("height_0.95_HPD", color = "#f8333c", size = 2, alpha = 0.5) +
	#geom_text(aes(label=node))+
	geom_text(aes(label=round(as.numeric(posterior), 2), x = branch), vjust=0, size = 2.5) +
#	geom_treescale(linesize = 1, offset = 2) +
        theme_tree2() + geom_rootedge(rootedge = 10) +
	scale_color_manual(values = ann_colors$ST, name = "ST") +
	scale_x_continuous(breaks = scales:::pretty_breaks(n = 5), minor_breaks = scales:::pretty_breaks(n = 15)) +
	theme(legend.position = "bottom", panel.grid.major = element_line(color = "black", size = .2),
	      panel.grid.minor = element_line(color = "grey", size = .2),
	      panel.grid.major.y = element_blank(),
              panel.grid.minor.y = element_blank()) +
	guides(colour = guide_legend(nrow = 1)) 

p1
ggsave(plot = p1, "", width = 8, height = 6)
