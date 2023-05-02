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
library(latex2exp)

# For the colours
colour <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6',
'#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3',
'#808000', '#ffd8b1', '#000075', '#808080', '#f0f0ff', '#000000')

ann_colors <- list(ST = c("Outgroup" = colour[3], "7" = colour[15], "8" = colour[10], "9" = colour[9], "10" = colour[6], "11" = colour[17], "12" = colour[8], "39" = colour[13], "40" = colour[2], "41" = colour[16], "42" = colour[7], "43" = colour[4], "71" = colour[11], "88" = colour[12], "102" = colour[19], "Geridu" = colour[5],"Ancient" = colour[1], "NF" = colour[20], "NIPH" = colour[22]))
clusterColours <- c("Western Mediterranean" = "#4D9DE0", "Fertile Crescent" = "#E15554", "Africa/America" = "#E1BC29", "Indo-Pacific" = "#3BB273", "Russia" = "#7768AE")
metaData <- read.delim("../AssemblyOnly/MetadataAll.tab", header =T) %>% as_tibble
dates <- read.delim("../BruceDatesUpdated.tab", header = F, col.names = c("Sample", "Date")) %>% as_tibble
metaData <- metaData %>% left_join(dates)

#### W. Med Phylogeny ####
# Basic Tree
tree <- read.beast(file = "AssemblyMLTrees/WMedSummary.nexus")
Tiplabels <- tree@phylo$tip.label

# Need to strip the dates off the ends of the tips for this tree
Tiplabels <- sapply(Tiplabels, function(x){
	tmp <- strsplit(x, split = "_") %>% unlist() %>% unname()
	tmp <- paste(tmp[1:(length(tmp)-1)], sep = "_", collapse = "_" )
	return(tmp)
	})
tree@phylo$tip.label <- Tiplabels

# Let's colour these based on their PhyloGroups 
STData <- read.delim("MLSTResults.txt")[,1:2]
rownames(STData) <- STData$Sample
STData["Brancorsini",] <- list(Sample = "Brancorsini", ST = "Ancient")
STData["Babortus-2308_1940","ST"] <-  "Outgroup"
STData["Reference",] <- list(Sample = "Reference", ST = "7")

STData <- STData[tree@phylo$tip.label,]
STData <- STData %>% mutate(ST = gsub("[^[:alnum:] ]","", ST))
STData <- STData %>% left_join(metaData %>% select(-ST)) %>% distinct() %>% 
	mutate(Type = ifelse(Sample == "Reference", "Reference", ifelse(Sample == "Brancorsini", "Ancient", "Modern"))) %>%
	mutate(Type = factor(Type, levels = c("Modern", "Ancient", "Reference")))

# Let's filter the ann_colors list so that only the STs present in the phylogeny are used
WMedColours <- ann_colors$ST[unique(STData$ST)]
WMedColours <- WMedColours[sort(names(WMedColours))]

# Getting the ranges in the phylogeny
tmp <- gsub("c\\(|\\)| ","",paste(tree@data$height_0.95_HPD,sep = ","))
tree@data$height_range_low <- as.numeric(gsub(".*,","",tmp))
tree@data$height_range_hi <- as.numeric(gsub(",.*","",tmp))

#  WMED
p1Wmed <- ggtree(tree, right = T, mrsd = "2017-01-01") %<+% STData +  #This is for the ST11 Phylo
	theme_tree2() +
	geom_rootedge(rootedge = 50) +
	geom_tiplab(aes(label = Country), align = T) + 
	annotate("rect", xmin = -8500, xmax = -6500, ymin = 2, ymax = 3.5, fill = "#3BB273", alpha = 0.5) + # Zeder 2008
	annotate("rect", xmin = -8500, xmax = -5800, ymin = 0.5, ymax = 2, fill = "#4D9DE0", alpha = 0.5) + # 
	annotate("rect", xmin = -7000, xmax = -5600, ymin = 3.5, ymax = 5, fill = "#7768AE", alpha = 0.5) + # 
	annotate("rect", xmin = -6400, xmax = -3700, ymin = 2, ymax = 3.5, fill = "#E1BC29", alpha = 0.5) +  # 
	annotate("text", x = -7500, y = 2.75, label = "Sheep\nDomestication") +
	annotate("text", x = -7150, y = 1.25, label = "Anatolia") +
	annotate("text", x = -6300, y = 4.25,label = "Mediterannean\nLittoral") + # Zeder 2008
	annotate("text", x = -5050, y = 2.75, label = "Central\nEurope") +
	geom_range("height_0.95_HPD", colour = "red", linewidth = 0.75, alpha = 0.75) +
	geom_tippoint(aes(colour = ST, shape = Type),size = 2) + 
	scale_color_manual(values = WMedColours) +
	guides(colour = guide_legend(nrow = 2, title.position = "top", title.hjust = 0.5), shape = guide_legend(nrow = 2, title.position = "top", title.hjust = 0.5, title = "Sample")) +
       	geom_nodelab(mapping = aes(label = ifelse(round(decimal_date(ymd("2017-01-01")) - height) < 1300,
							 paste0("[",paste(round(decimal_date(ymd("2017-01-01")) - height_range_low),round(decimal_date(ymd("2017-01-01")) - height_range_hi), sep = ","), "]"), NA)), # This is for the ST11 Phylo
			    geom = "text", nudge_y = 0.4, nudge_x = -250) +
	xlab("Year") +
	scale_x_continuous(breaks = scales:::pretty_breaks(n = 10)) +
	theme(panel.grid.major.x = element_line(color = "grey50", size = 0.2),
	      panel.grid.minor.x = element_line(color = "grey90", size = 0.2),
	      legend.position = "bottom") 

treeData <- tree@data[,c(10,1:3)]

translatedDates <- treeData %>% mutate(height = as.numeric(height),ymd =decimal_date(ymd("2017-01-01")),  heightNew = round(decimal_date(ymd("2017-01-01")) - height))

translatedDates <- treeData %>% mutate(height = as.numeric(height),height = round(decimal_date(ymd("2017-01-01")) - height),height_median = round(decimal_date(ymd("2017-01-01")) - height_median),height_0.95_HPD = gsub("c\\(|\\)| ","",paste(height_0.95_HPD,sep = ","))) %>%
       	filter(height_median < 1500) %>% separate(height_0.95_HPD, into = c("HeightCIHi", "HeightCILo"), sep = ",") %>% mutate(HeightCILo = round(decimal_date(ymd("2017-01-01")) - as.numeric(HeightCILo))) %>%
	mutate(HeightCIHi = round(decimal_date(ymd("2017-01-01")) - as.numeric(HeightCIHi)))

ggsave(p1Wmed, file = "WmedPhylogenyML.pdf", width = 9, height = 6)

#### W. Med Phylogeny BEAST ####
# Basic Tree
tree <- read.beast(file = "Beast/SummarizedTrees/WMedStrict.nexus")
Tiplabels <- tree@phylo$tip.label

# Need to strip the dates off the ends of the tips for this tree
#Tiplabels <- sapply(Tiplabels, function(x){
#	tmp <- strsplit(x, split = "_") %>% unlist() %>% unname()
#	tmp <- paste(tmp[1:(length(tmp)-1)], sep = "_", collapse = "_" )
#	return(tmp)
#	})
#tree@phylo$tip.label <- Tiplabels

# Let's colour these based on their PhyloGroups 
STData <- read.delim("MLSTResults.txt")[,1:2]
rownames(STData) <- STData$Sample
STData["Brancorsini",] <- list(Sample = "Brancorsini", ST = "Ancient")
STData["Babortus-2308_1940","ST"] <-  "Outgroup"
STData["Reference",] <- list(Sample = "Reference", ST = "7")

STData <- STData[tree@phylo$tip.label,]
STData <- STData %>% mutate(ST = gsub("[^[:alnum:] ]","", ST))
STData <- STData %>% left_join(metaData %>% select(-ST)) %>% distinct() %>% 
	mutate(Type = ifelse(Sample == "Reference", "Reference", ifelse(Sample == "Brancorsini", "Ancient", "Modern"))) %>%
	mutate(Type = factor(Type, levels = c("Modern", "Ancient", "Reference")))

# Let's filter the ann_colors list so that only the STs present in the phylogeny are used
WMedColours <- ann_colors$ST[unique(STData$ST)]
WMedColours <- WMedColours[sort(names(WMedColours))]

# Getting the ranges in the phylogeny
tmp <- gsub("c\\(|\\)| ","",paste(tree@data$height_0.95_HPD,sep = ","))
tree@data$height_range_low <- as.numeric(gsub(".*,","",tmp))
tree@data$height_range_hi <- as.numeric(gsub(",.*","",tmp))

#  WMED
p1Wmed <- ggtree(tree, right = T, mrsd = "2017-01-01") %<+% STData +  #This is for the ST11 Phylo
	theme_tree2() +
	geom_rootedge(rootedge = 50) +
	geom_tiplab(aes(label = Country), align = T) + 
	annotate("rect", xmin = -8500, xmax = -6500, ymin = 2, ymax = 3.5, fill = "#3BB273", alpha = 0.5) + # Zeder 2008
	annotate("rect", xmin = -8500, xmax = -5800, ymin = 0.5, ymax = 2, fill = "#4D9DE0", alpha = 0.5) + # 
	annotate("rect", xmin = -7000, xmax = -5600, ymin = 3.5, ymax = 5, fill = "#7768AE", alpha = 0.5) + # 
	annotate("rect", xmin = -6400, xmax = -3700, ymin = 2, ymax = 3.5, fill = "#E1BC29", alpha = 0.5) +  # 
	annotate("text", x = -7500, y = 2.75, label = "Sheep\nDomestication") +
	annotate("text", x = -7150, y = 1.25, label = "Anatolia") +
	annotate("text", x = -6300, y = 4.25,label = "Mediterannean\nLittoral") + # Zeder 2008
	annotate("text", x = -5050, y = 2.75, label = "Central\nEurope") +
	geom_range("height_0.95_HPD", colour = "red", linewidth = 0.75, alpha = 0.75) +
	geom_tippoint(aes(colour = ST, shape = Type),size = 2) + 
	scale_color_manual(values = WMedColours) +
	guides(colour = guide_legend(nrow = 2, title.position = "top", title.hjust = 0.5), shape = guide_legend(nrow = 2, title.position = "top", title.hjust = 0.5, title = "Sample")) +
       	geom_nodelab(mapping = aes(label = ifelse(round(decimal_date(ymd("2017-01-01")) - height) < 1300,
							 paste0("[",paste(round(decimal_date(ymd("2017-01-01")) - height_range_low),round(decimal_date(ymd("2017-01-01")) - height_range_hi), sep = ","), "]"), NA)), # This is for the ST11 Phylo
			    geom = "text", nudge_y = 0.4, nudge_x = -250) +
	xlab("Year") +
	scale_x_continuous(breaks = scales:::pretty_breaks(n = 10)) +
	theme(panel.grid.major.x = element_line(color = "grey50", size = 0.2),
	      panel.grid.minor.x = element_line(color = "grey90", size = 0.2),
	      legend.position = "bottom") 

treeData <- tree@data[,c(10,1:3)]

translatedDates <- treeData %>% mutate(height = as.numeric(height),ymd =decimal_date(ymd("2017-01-01")),  heightNew = round(decimal_date(ymd("2017-01-01")) - height))

translatedDates <- treeData %>% mutate(height = as.numeric(height),height = round(decimal_date(ymd("2017-01-01")) - height),height_median = round(decimal_date(ymd("2017-01-01")) - height_median),height_0.95_HPD = gsub("c\\(|\\)| ","",paste(height_0.95_HPD,sep = ","))) %>%
       	filter(height_median < 1500) %>% separate(height_0.95_HPD, into = c("HeightCIHi", "HeightCILo"), sep = ",") %>% mutate(HeightCILo = round(decimal_date(ymd("2017-01-01")) - as.numeric(HeightCILo))) %>%
	mutate(HeightCIHi = round(decimal_date(ymd("2017-01-01")) - as.numeric(HeightCIHi)))

ggsave(p1Wmed, file = "WmedPhylogenyBEAST.pdf", width = 9, height = 6)

##### WMed with Geridu #####

tree <- read.newick(file = "NewMLTrees/GeriduWMed.treefile")
tree <- phytools::reroot(tree, interactive = T)
Tiplabels <- tree$tip.label

# Need to strip the dates off the ends of the tips for this tree
Tiplabels <- sapply(Tiplabels, function(x){
	tmp <- strsplit(x, split = "_") %>% unlist() %>% unname()
	tmp <- paste(tmp[1:(length(tmp)-1)], sep = "_", collapse = "_" )
	return(tmp)
	})

ggtree(tree, right = T) +
	geom_tiplab() +
	geom_treescale(linesize = 1, offset = 2) + theme_tree()# + geom_rootedge(rootedge = 0.00005) 
ggsave("GeriduPhylogeny.pdf", width = 9, height = 6)

#### BEAST Global Phylogeny ####
# Basic Tree
tree <- read.beast(file = "Beast/SummarizedTrees/GlobalStrict.nexus")
Tiplabels <- tree@phylo$tip.label

# Need to strip the dates off the ends of the tips for this tree
#Tiplabels <- sapply(Tiplabels, function(x){
#	tmp <- strsplit(x, split = "_") %>% unlist() %>% unname()
#	tmp <- paste(tmp[1:(length(tmp)-1)], sep = "_", collapse = "_" )
#	return(tmp)
#	})
#tree@phylo$tip.label <- Tiplabels

# Let's colour these based on their PhyloGroups 
STData <- read.delim("MLSTResults.txt")[,1:2]
rownames(STData) <- STData$Sample
STData["Brancorsini",] <- list(Sample = "Brancorsini", ST = "Ancient")
STData["Babortus-2308_1940","ST"] <-  "Outgroup"
STData["Reference",] <- list(Sample = "Reference", ST = "7")
STData <- STData %>% mutate(ST = gsub("[^[:alnum:] ]","", ST)) %>% left_join(metaData %>% select(-ST)) %>% distinct() 

# Now to label the NIPH Strains
STData <- STData %>% mutate(ST = replace(ST, grepl("^NIPH-*|NI_2007", Sample) & grepl("8",ST), "NIPH")) %>%
	mutate(Type = ifelse(Sample == "Reference", "Reference", ifelse(Sample == "Brancorsini", "Ancient", "Modern"))) %>%
	mutate(Type = factor(Type, levels = c("Modern", "Ancient", "Reference"))) %>% as.data.frame() %>%
	filter(!grepl("BY38|Brucellamelitensis-043", Sample))
rownames(STData) <- STData$Sample

STData <- STData[tree@phylo$tip.label,]

# Let's filter the ann_colors list so that only the STs present in the phylogeny are used
GlobalColours <- ann_colors$ST[unique(STData$ST)]
GlobalColours <- GlobalColours[mixedsort(names(GlobalColours))]

tmp <- gsub("c\\(|\\)| ","",paste(tree@data$height_0.95_HPD,sep = ","))
tree@data$height_range_low <- as.numeric(gsub(".*,","",tmp))
tree@data$height_range_hi <- as.numeric(gsub(",.*","",tmp))

# Global
p1Global <- ggtree(tree, right = T, mrsd = "2018-06-01") %<+% STData +  #This is for the ST11 Phylo
	theme_tree2() +
	#geom_text2(aes(subset=!isTip, label = node)) +
	geom_rootedge(rootedge = 50) +
	annotate("rect", xmin = -8500, xmax = -6500, ymin = 4.5, ymax = 8.5, fill = "#3BB273", alpha = 0.5) + # Zeder 2008
	annotate("rect", xmin = -8500, xmax = -5800, ymin = 0.5, ymax = 4.5, fill = "#4D9DE0", alpha = 0.5) + # 
	annotate("rect", xmin = -7000, xmax = -5600, ymin = 8.5, ymax = 12.5, fill = "#7768AE", alpha = 0.5) + # 
	annotate("rect", xmin = -6400, xmax = -3700, ymin = 4.5, ymax = 8.5, fill = "#E1BC29", alpha = 0.5) +  # 
	annotate("text", x = -7500, y = 6.5, label = "Sheep\nDomestication") +
	annotate("text", x = -7150, y = 2.5, label = "Anatolia") +
	annotate("text", x = -6300, y = 10.5,label = "Mediterannean\nLittoral") + # Zeder 2008
	annotate("text", x = -5050, y = 6.5, label = "Central\nEurope") +
	geom_range("height_0.95_HPD", colour = "red", linewidth = 0.75, alpha = 0.75) +
	geom_tippoint(aes(colour = ST, shape = Type),size = 2) + 
	scale_color_manual(values = GlobalColours) +
	guides(colour = guide_legend(nrow = 2, title.position = "top", title.hjust = 0.5), shape = guide_legend(nrow = 2, title.position = "top", title.hjust = 0.5, title = "Sample")) +
       	geom_nodelab(mapping = aes(label = ifelse(round(decimal_date(ymd("2018-06-01")) - height) < 1300,
							 paste0("[",paste(round(decimal_date(ymd("2018-06-01")) - height_range_low),round(decimal_date(ymd("2017-01-01")) - height_range_hi), sep = ","), "]"), NA)), # This is for the ST11 Phylo
			    geom = "text", nudge_y = 0.4, nudge_x = -250) +
	xlab("Year") +
	scale_x_continuous(breaks = scales:::pretty_breaks(n = 10)) +
	theme(panel.grid.major.x = element_line(color = "grey50", size = 0.2),
	      panel.grid.minor.x = element_line(color = "grey90", size = 0.2),
	      legend.position = "bottom") 

# Going to try replacing the labels with the Phylogroups
p1Global <- p1Global + 
	geom_cladelabel(align = T, node = 122, label = "Western\nMediterranean", colour = ann_colors$ST["11"]) +
	geom_cladelabel(align = T, node = 146, label = "Eastern\nMediterranean", colour = ann_colors$ST["8"]) +
	geom_cladelabel(align = T, node = 228, label = "African", colour = ann_colors$ST["12"]) +
	geom_cladelabel(align = T, node = 222, label = "American", colour = ann_colors$ST["7"]) 

ggsave(p1Global, file = "GlobalPhylogenyBeast.pdf", width = 9, height = 6)


treeData <- tree@data[,c(14,1:3)]

translatedDates <- treeData %>% mutate(height = as.numeric(height),ymd =decimal_date(ymd("2018-06-01")),  heightNew = round(decimal_date(ymd("2018-06-01")) - height))
translatedDates <- treeData %>% mutate(height = as.numeric(height),height = round(decimal_date(ymd("2018-06-01")) - height),height_median = round(decimal_date(ymd("2018-06-01")) - height_median),height_0.95_HPD = gsub("c\\(|\\)| ","",paste(height_0.95_HPD,sep = ","))) %>%
       	filter(height_median < 1300) %>% separate(height_0.95_HPD, into = c("HeightCIHi", "HeightCILo"), sep = ",") %>% mutate(HeightCILo = round(decimal_date(ymd("2018-06-01")) - as.numeric(HeightCILo))) %>%
	mutate(HeightCIHi = round(decimal_date(ymd("2017-01-01")) - as.numeric(HeightCIHi)))

##### Global Phylogeny #####
# Basic Tree
tree <- read.newick(file = "AssemblyMLTrees/GlobalReview.treefile")
tree <- phytools::reroot(tree, interactive = T)
Tiplabels <- tree$tip.label

# Let's colour these based on their PhyloGroups 
STData <- read.delim("MLSTResults.txt")[,1:2]
rownames(STData) <- STData$Sample
STData["Brancorsini",] <- list(Sample = "Brancorsini", ST = "Ancient")
STData["Babortus-2308_1940","ST"] <-  "Outgroup"
STData["Reference",] <- list(Sample = "Reference", ST = "7")
STData <- STData %>% mutate(ST = gsub("[^[:alnum:] ]","", ST))

# Now to label the NIPH Strains
STData <- STData %>% mutate(ST = replace(ST, grepl("^NIPH-*|NI_2007", Sample) & grepl("8",ST), "NIPH")) %>%
	mutate(Type = ifelse(Sample == "Reference", "Reference", ifelse(Sample == "Brancorsini", "Ancient", "Modern"))) %>%
	mutate(Type = factor(Type, levels = c("Modern", "Ancient", "Reference")))

STData <- STData[tree$tip.label,]

# Let's filter the ann_colors list so that only the STs present in the phylogeny are used
GlobalColours <- ann_colors$ST[unique(STData$ST)]
GlobalColours <- GlobalColours[mixedsort(names(GlobalColours))]

# Going to try replacing the labels with the Phylogroups
p1Global <- ggtree(tree, right = T) %<+% STData +
	#geom_text2(aes(subset=!isTip, label = node)) +
       	#geom_nodepoint(size = 1.5, colour = ifelse(as.numeric(tree$node.label) >= 90, "black", ifelse(as.numeric(tree$node.label) >= 50, "grey",NA)),
	#	     shape = "square") +
	geom_tippoint(mapping = aes(colour = ST, shape = Type), size = 1.5) +
	#geom_tiplab(align = T) + 
	geom_treescale(linesize = 1, offset = 2) + theme_tree() + geom_rootedge(rootedge = 0.00005) +
	scale_color_manual(values = GlobalColours, name = "Sequence Type") +
	theme(legend.position = "bottom") + guides(colour = guide_legend(nrow = 2, title.position = "top", title.hjust = 0.5), shape = guide_legend(nrow = 2, title.position = "top", title.hjust = 0.5, title = "Sample"))

p1Global <- p1Global + 
	geom_cladelabel(align = T, node = 121, label = "Western\nMediterranean", colour = ann_colors$ST["11"]) +
	geom_cladelabel(align = T, node = 145, label = "Eastern\nMediterranean", colour = ann_colors$ST["8"]) +
	geom_cladelabel(align = T, node = 228, label = "African", colour = ann_colors$ST["12"]) +
	geom_cladelabel(align = T, node = 221, label = "American", colour = ann_colors$ST["7"]) 

ggsave(p1Global, file = "GlobalPhylogeny.pdf", width = 9, height = 6)


# Now to see if the Accessory Clustering reveals anything interesting
clusterInfo <- read.delim("../PanGenomeAnalysisFixedMapping/GenomeAccessClustered.tab") %>%
	mutate(clusters = ifelse(clusters == 1, "Western Mediterranean",
				 ifelse(clusters == 2, "Fertile Crescent",
				       	ifelse(clusters == 3, "Africa/America", ifelse(clusters == 4, "Indo-Pacific", "Russia")))))
clusterInfo[which(clusterInfo$Genome == "Geridu"),"Genome"] <- "Nodule1_S1"
rownames(clusterInfo) <- clusterInfo$Genome
clusterInfo <- clusterInfo[tree$tip.label,]

# Seeing if there's anything going on with the clusters
p1Cluster <- ggtree(tree, right = T) %<+% clusterInfo +
	#geom_text2(aes(subset=!isTip, label = node)) +
       	#geom_nodepoint(size = 1.5, colour = ifelse(as.numeric(tree$node.label) >= 90, "black", ifelse(as.numeric(tree$node.label) >= 50, "grey",NA)),
	#	     shape = "square") +
	geom_tippoint(mapping = aes(colour = as.factor(clusters)), size = 1.5) +
	#geom_tiplab(align = T) + 
	geom_treescale(linesize = 1, offset = 2) + theme_tree() + geom_rootedge(rootedge = 0.00005) +
	scale_color_manual(values = clusterColours, name = "Accessory Genome Clustering") +
	theme(legend.position = "bottom") + guides(colour = guide_legend(nrow = 2))

ggarrange(p1Global, p1Wmed, common.legend = T, legend = "bottom", nrow = 2, labels = "auto", align = "hv")
ggsave(file = "MLPhylogenies.pdf", width = 12, height = 9)

ggsave(p1Global, file = "GlobalML.pdf", width = 9, height = 6)
ggsave(p1Cluster, file = "ClusterML.pdf", width = 9, height = 6)
ggsave(p1Wmed, file = "WMedML.pdf", width = 9, height = 6)

##### E. Med Phylogeny #####
tree <- read.newick(file = "NewMLTrees/EMed.treefile")
tree <- phytools::reroot(tree, interactive = T)
Tiplabels <- tree$tip.label

# Let's colour these based on their PhyloGroups 
STData <- read.delim("MLSTResults.txt")[,1:2]
rownames(STData) <- STData$Sample
STData["Brancorsini",] <- list(Sample = "Brancorsini", ST = "Brancorsini")
STData["Reference",] <- list(Sample = "Reference", ST = "7")
STData["Nodule1_S1",] <- list(Sample = "Nodule1_S1", ST = "Geridu")

STData <- STData[tree$tip.label,]
STData <- STData %>% mutate(ST = gsub("[^[:alnum:] ]","", ST))

clusterInfo <- read.delim("../PanGenomeAnalysisFixedMapping/GenomeAccessClustered.tab") %>%
	mutate(clusters = ifelse(clusters == 1, "Western Mediterranean",
				 ifelse(clusters == 2, "Fertile Crescent",
				       	ifelse(clusters == 3, "Africa/America", ifelse(clusters == 4, "Indo-Pacific", "Russia")))))

clusterInfo[which(clusterInfo$Genome == "Geridu"),"Genome"] <- "Nodule1_S1"
rownames(clusterInfo) <- clusterInfo$Genome
clusterInfo <- clusterInfo[tree$tip.label,] %>% left_join(metaData, by = c("Genome" = "Sample"))

# Seeing if there's anything going on with the clusters
p1Cluster <- ggtree(tree, right = T) %<+% clusterInfo +
	geom_tippoint(mapping = aes(colour = as.factor(clusters)), size = 1.5) +
	geom_tiplab(align = T, mapping = aes(label = Country, colour = clusters), show.legend = F) +
	#geom_tiplab(align = T) + 
	geom_treescale(linesize = 1, offset = 2) + theme_tree() + geom_rootedge(rootedge = 0.00005) +
	scale_color_manual(values = clusterColours, name = "Accessory Genome Clustering") +
	theme(legend.position = "bottom") + guides(colour = guide_legend(nrow = 2))

ggsave(p1Cluster, file = "EmedClusterML.pdf", height = 12, width = 12)
