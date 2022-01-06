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
library(lubridate)

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
	theme_tree2() +
	geom_rootedge(rootedge = 50) +
	geom_range("height_0.95_HPD", colour = "red", size = 0.75, alpha = 0.75) +
	geom_tippoint(aes(colour = ST)) +
	scale_color_manual(values = ann_colors$ST, name = "Sequence Type") +
	new_scale_color() +
       	geom_nodepoint(aes(color = ifelse(posterior < 0.5, NA,
						 ifelse(posterior >= 0.5 & posterior < 0.9, "Fifty", "Ninety"))),
			      shape = "square", show.legend = F) +
	scale_color_manual(values = c("NA" = NA, "Fifty" = "grey", "Ninety" = "black")) +
       	geom_nodelab(size = 2.5,mapping = aes(label = round(2017.0014 - height_median)),
			    geom = "label", nudge_y = 0.4, nudge_x = -50) +
	xlab("Year") +
	scale_x_continuous(breaks = scales:::pretty_breaks()) +
	theme(panel.grid.major.x = element_line(color = "grey10", size = 0.2),
	      panel.grid.minor.x = element_line(color = "grey80", size = 0.2),
	      legend.position = "bottom")# + guides(colour = guide_legend(nrow = 1))

ggsave(plot = p1, "FirstBeastST11Trimmed.pdf", width = 12, height = 9)

decimal_date(ymd("1394-08-14"))
