library(dplyr)
library(tidyr)
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

ann_colors <- list(ST = c("5" = colour[10], "7" = colour[15], "8" = colour[3], "9" = colour[9], "10" = colour[6], "11" = colour[17], "12" = colour[8], "39" = colour[13], "40" = colour[2], "41" = colour[16], "42" = colour[7], "43" = colour[4], "71" = colour[11], "88" = colour[12], "102" = colour[19], "Brancorsini" = colour[1], "NF" = colour[20], "Reference" = colour[22], "Geridu" = colour[5]))
#		Clusters = c("1" = colour[4],"2" = colour[2], "3" = colour[5], "4" = colour[6], "5" = colour[1], "0" = colour[22]))

LocationID <- function(Sample){
	ifelse(grepl("\\dTE\\d|\\dCB\\d|\\dIS\\d|Kay|JessSample", Sample), "Italy",
	       ifelse(grepl("BwIM-MAR-25", Sample), "Morocco",
		      ifelse(grepl("CT-US", Sample), "US",
			     ifelse(grepl("UK", Sample), "UK",
				    ifelse(grepl("BRUC101", Sample), "Egypt",
					   ifelse(grepl("00-2529-3", Sample),"France", "Unknown"))))))
}

# Basic Tree
tree <- read.beast(file = "BeastResults/InvariantWhole/WholeInvariant.nexus")
Tiplabels <- tree@phylo$tip.label

####################################################
# Let's colour these based on their Sequence Type ##
####################################################
STData <- read.delim("../MLSTResults.txt")[,1:2]
rownames(STData) <- STData$Sample
STData["JessSample",] <- list(Sample = "JessSample", ST = "Brancorsini")
STData["Reference",] <- list(Sample = "Reference", ST = "Reference")
STData["KayBMel",] <- list(Sample = "KayBMel", ST = "Geridu")

STData<- STData[Tiplabels,]
STData <- STData %>% mutate(ST = gsub("[^[:alnum:] ]","", ST))

# Let's filter the ann_colors list so that only the STs present in the phylogeny are used
ann_colors$ST <- ann_colors$ST[mixedsort(unique(STData$ST))]

# Let's see if we can include some geographic information
STData <- STData %>% mutate(Country = LocationID(Sample))

# Reading the Clustering data from the P/A Analysis
#wholeCluster <- read.delim("../PanGenomeAnalysis/ItalyClustering.tab") %>% as_tibble()
#colnames(wholeCluster)[5] <- "Sample"
#wholeCluster$Sample[which(wholeCluster$Sample %in% c("Brancorsini", "Geridu"))] <- c("JessSample", "KayBMel")
#STData <- STData %>% left_join(wholeCluster %>% select(Sample, clusters)) %>% mutate(clusters = replace(clusters, is.na(clusters), 0))

############Cluster Plotting#################
##pClust <- ggtree(tree, right = T, mrsd = "2018-05-01") %<+% STData +
#pClust <- ggtree(tree, right = T, mrsd = "2017-01-01") %<+% STData + #This is for the ST11 Phylo
#	theme_tree2() +
#	geom_rootedge(rootedge = 50) +
#	geom_range("height_0.95_HPD", colour = "red", size = 0.75, alpha = 0.75) +
#	#geom_tippoint(aes(colour = ST), size = 2) +
#	geom_tippoint(aes(colour = as.factor(clusters)), size = 2) + #only if including clustering results
#	#geom_tippoint(aes(colour = ST, shape = Country), size = 2) + #only for ST11
#	scale_color_manual(values = ann_colors$Clusters) + guides(colour = guide_legend(nrow = 2, title ="Clusters")) +
#	#scale_shape_manual(values = c(16,15,17,18,10,14,8)) + # only for ST11
#	#scale_shape_manual(values = c(16:18,6:10)) +
#	new_scale_color() +
#       	geom_nodepoint(aes(color = ifelse(posterior < 0.5, NA,
#						 ifelse(posterior >= 0.5 & posterior < 0.9, "Fifty", "Ninety"))),
#			      shape = "square", show.legend = F) +
#	scale_color_manual(values = c("NA" = NA, "Fifty" = "grey", "Ninety" = "black")) +
#       	#geom_nodelab(size = 2.5,mapping = aes(label = round(2018.3287671232877 - height_median)),
#       	#geom_nodelab(size = 2.5,mapping = aes(label = round(decimal_date(ymd("2017-01-01")) - height_median)), # This is for the ST11 Phylo
#	#		    geom = "label", nudge_y = 0.4, nudge_x = -50) +
#	xlab("Year") +
#	scale_x_continuous(breaks = scales:::pretty_breaks()) +
#	theme(panel.grid.major.x = element_line(color = "grey10", size = 0.2),
#	      panel.grid.minor.x = element_line(color = "grey80", size = 0.2),
#	      legend.position = "bottom")
#
#pClust
#ggsave("Beast/ST11Clustered.pdf", width = 9, height = 6)

p1 <- ggtree(tree, right = T, mrsd = "2018-05-01") %<+% STData +
#p1 <- ggtree(tree, right = T, mrsd = "2017-01-01") %<+% STData + #This is for the ST11 Phylo
	theme_tree2() +
	geom_rootedge(rootedge = 50) +
	geom_range("height_0.95_HPD", colour = "red", size = 0.75, alpha = 0.75) +
	geom_tippoint(aes(colour = ST), size = 2) +
	#geom_tippoint(aes(colour = ST, shape = Country), size = 2) + #only for ST11
	scale_color_manual(values = ann_colors$ST) + guides(colour = guide_legend(nrow = 2, title ="Sequence Type")) +
	#scale_shape_manual(values = c(16,15,17,18,10,14,8)) + # only for ST11
	#scale_shape_manual(values = c(16:18,6:10)) +
	new_scale_color() +
       	geom_nodepoint(aes(color = ifelse(posterior < 0.5, NA,
						 ifelse(posterior >= 0.5 & posterior < 0.9, "Fifty", "Ninety"))),
			      shape = "square", show.legend = F) +
	scale_color_manual(values = c("NA" = NA, "Fifty" = "grey", "Ninety" = "black")) +
       	geom_nodelab(size = 2.5,mapping = aes(label = round(2018.3287671232877 - height_median)),
       	#geom_nodelab(size = 2.5,mapping = aes(label = round(decimal_date(ymd("2017-01-01")) - height_median)), # This is for the ST11 Phylo
			    geom = "label", nudge_y = 0.4, nudge_x = -50) +
	xlab("Year") +
	scale_x_continuous(breaks = scales:::pretty_breaks()) +
	theme(panel.grid.major.x = element_line(color = "grey10", size = 0.2),
	      panel.grid.minor.x = element_line(color = "grey80", size = 0.2),
	      legend.position = "bottom")

p1 

ggsave("~/Documents/University/LabMeetings/2022/Mar25/Figures/WholePhylo.pdf", width = 9, height = 6)

tree@data %>% mutate(height_median = round(decimal_date(ymd("2017-01-01")) - height_median),heightn = round(decimal_date(ymd("2017-01-01")) - height),height_0.95_HPD = gsub("c\\(|\\)| ","",paste(height_0.95_HPD,sep = ","))) %>%
       	filter(height_median < 1500) %>% separate(height_0.95_HPD, into = c("HeightCIHi", "HeightCILo"), sep = ",") %>% mutate(HeightCILo = round(decimal_date(ymd("2017-01-01")) - as.numeric(HeightCILo))) %>%
	mutate(HeightCIHi = round(decimal_date(ymd("2017-01-01")) - as.numeric(HeightCIHi)))

#ggsave(plot = p1, "ST11Invariant.pdf", width = 9, height = 6)
#######################
### Now for Tempest ###
#######################

# Now to get the Tempest Plot
tempest <- as_tibble(read.delim("BeastResults/InvariantST11/ST11KayIncludedInvariantTempestSNPS.tab"))# %>% filter(date > 0)
# Getting the Phylogroups plotted

STData <- STData[match(tempest$tip, STData$Sample),] # Getting them in the Right Order
tempest$ST <- STData$ST
#tempest$Country <- STData$Country

tempestPlot <- tempest %>% ggplot(aes(x = date, y  = distance)) + geom_smooth(method = "lm", show.legend = F, colour = "black") +
	geom_point(alpha = 0.75, aes(colour = ST), show.legend = T) + theme_bw() + ylab("Root to Tip Divergence") + xlab("Year") +
	#scale_shape_manual(values = c(16,15,17,18,10,14,8)) +
	scale_color_manual(values = ann_colors$ST, name = "Sequence Type")
tempestPlot

model <- lm(distance ~date, data = tempest)
summary(model)
r2 <- round(summary(model)$adj.r.squared,3)
tmp <- summary(model)$fstatistic
pval <- round(pf(tmp[1],tmp[2],tmp[3], lower.tail = F),3)

tempestPlot <- tempestPlot + annotate(geom = "text", x = 1700, y = 0.0026,label = bquote(R[adj]^2 == .(r2))) +
	annotate(geom = "text", x = 1700, y = 0.00255, label = bquote(P == .(pval))) +
	scale_y_continuous(breaks = scales:::pretty_breaks(n = 10)) + theme(legend.position = "bottom") + #ggtitle("ST11 Clade") +
	guides(colour = guide_legend(nrow = 1))
ggsave("~/Documents/University/LabMeetings/2022/Mar25/Figures/ST11Tempest.pdf", width = 6, height = 4)

ggarrange(p1,tempestPlot, ncol = 1, common.legend = T, legend = "bottom", labels = "AUTO", heights = c(2,1))
ggsave(file = "WholePhylogeneticsSNPsWithClusters.pdf", width = 9, height = 12)

