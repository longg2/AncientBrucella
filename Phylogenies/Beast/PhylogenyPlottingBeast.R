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
library(xlsx)

# For the colours
colour <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6',
'#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3',
'#808000', '#ffd8b1', '#000075', '#808080', '#f0f0ff', '#000000')

ann_colors <- list(ST = c("5" = colour[10], "7" = colour[15], "8" = colour[3], "9" = colour[9], "10" = colour[6], "11" = colour[17], "12" = colour[8], "39" = colour[13], "40" = colour[2], "41" = colour[16], "42" = colour[7], "43" = colour[4], "71" = colour[11], "88" = colour[12], "102" = colour[19], "Ancient" = colour[1], "NF" = colour[20], "Reference" = colour[22], "Kay" = colour[5]))

LocationID <- function(Sample){
	ifelse(grepl("\\dTE\\d|\\dCB\\d|\\dIS\\d|Kay|JessSample", Sample), "Italy",
	       ifelse(grepl("BwIM-MAR-25", Sample), "Morocco",
		      ifelse(grepl("CT-US", Sample), "US",
			     ifelse(grepl("UK", Sample), "UK",
				    ifelse(grepl("BRUC101", Sample), "Egypt",
					   ifelse(grepl("00-2529-3", Sample),"France", "Unknown"))))))
}

# Basic Tree
tree <- read.beast(file = "BeastResults/ST11Strict/test.nexus")
Tiplabels <- tree@phylo$tip.label

####################################################
# Let's colour these based on their Sequence Type ##
####################################################
STData <- read.delim("../MLSTResults.txt")[,1:2]
rownames(STData) <- STData$Sample
STData["JessSample",] <- list(Sample = "JessSample", ST = "Ancient")
STData["Reference",] <- list(Sample = "Reference", ST = "Reference")
STData["KayBMel",] <- list(Sample = "KayBMel", ST = "Kay")

STData<- STData[Tiplabels,]
STData <- STData %>% mutate(ST = gsub("[^[:alnum:] ]","", ST))

# Let's filter the ann_colors list so that only the STs present in the phylogeny are used
ann_colors$ST <- ann_colors$ST[unique(STData$ST)]

# Let's see if we can include some geographic information
STData <- STData %>% mutate(Country = LocationID(Sample))

#p1 <- ggtree(tree, right = T, mrsd = "2018-05-01") %<+% STData +
p1 <- ggtree(tree, right = T, mrsd = "2017-01-01") %<+% STData + #This is for the ST11 Phylo
	theme_tree2() +
	geom_rootedge(rootedge = 50) +
	geom_range("height_0.95_HPD", colour = "red", size = 0.75, alpha = 0.75) +
	geom_tippoint(aes(colour = ST, shape = Country), size = 3) +
	scale_color_manual(values = ann_colors$ST, name = "Sequence Type") +
	scale_shape_manual(values = c(16:18,6:10)) +
	new_scale_color() +
       	geom_nodepoint(aes(color = ifelse(posterior < 0.5, NA,
						 ifelse(posterior >= 0.5 & posterior < 0.9, "Fifty", "Ninety"))),
			      shape = "square", show.legend = F) +
	scale_color_manual(values = c("NA" = NA, "Fifty" = "grey", "Ninety" = "black")) +
       	#geom_nodelab(size = 2.5,mapping = aes(label = round(2018.3287671232877 - height_median)),
       	geom_nodelab(size = 2.5,mapping = aes(label = round(2017.0014 - height_median)), # This is for the ST11 Phylo
			    geom = "label", nudge_y = 0.4, nudge_x = -50) +
	xlab("Year") +
	scale_x_continuous(breaks = scales:::pretty_breaks()) +
	theme(panel.grid.major.x = element_line(color = "grey10", size = 0.2),
	      panel.grid.minor.x = element_line(color = "grey80", size = 0.2),
	      legend.position = "bottom")# + guides(colour = guide_legend(nrow = 1))

	p1

ggsave(plot = p1, "ST11BeastWithKay.pdf", width = 12, height = 9)

decimal_date(ymd("1394-08-14"))
