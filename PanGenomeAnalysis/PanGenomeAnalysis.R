library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(ggvenn)
library(scales)
library(ggExtra)
library(reshape2)
library(cluster)
library(parallel)
library(pbapply)
library(purrr)
library(pheatmap)
library(dendextend)
library(FactoMineR)
library(factoextra)
library(gtools)

## Functions
CorrectingGC <- function(depthData){ # From the iREP Paper
	GCModel <- lm(Mean ~ GCContent, data = depthData)
	resids <- abs(summary(GCModel)$residuals) # Trying to find the top 1% largest residuals. Don't care about location
	residFilt <- resids < quantile(resids,0.99)
	GCModelFilt <- lm(Mean ~ GCContent, data = depthData[residFilt,])
	adjr2 <- summary(GCModel)$adj.r.squared
	beforeGC <- depthData[residFilt,] %>% ggplot(aes(x = GCContent, y = Mean)) +
		geom_point() + geom_density_2d(colour = "grey") + theme_bw() + geom_smooth(method = "lm") +
	       	ggtitle(bquote("Before GC Correction"~R[adj]^2 == .(round(adjr2, 3))), subtitle = "Top 1% Residuals Filtered") + 
	       	ylab("Mean Coverage") + xlab("GC Content")
	if(adjr2 > 0){ # We make the correction
		covAverage <- mean(depthData$Mean)
		correction <- covAverage - predict(GCModel, data = depthData$GCContent)
		depthData <- depthData %>% mutate(GCCorrected = Mean + correction)
	}else{
		depthData <- depthData %>% mutate(GCCorrected = Mean)
	}
	NewGCModel <- lm(GCCorrected ~ GCContent, data = depthData[residFilt,])
	adjr2 <- summary(NewGCModel)$adj.r.squared
	afterGC <- depthData[residFilt,] %>% ggplot(aes(x = GCContent, y = GCCorrected)) +
		geom_point() + geom_density_2d(colour = "grey") + theme_bw() + geom_smooth(method = "lm") +
	       	ggtitle(bquote("After GC Correction"~R[adj]^2 == .(round(adjr2, 3))), subtitle = "Top 1% Residuals Filtered") +
	       	ylab("Mean Coverage") + xlab("GC Content")
	GCPlots <- ggarrange(beforeGC, afterGC, ncol = 2, align = "hv")
	return(list("Plot" = GCPlots, UpdatedDepths = depthData))

}
############

colour <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6',
'#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3',
'#808000', '#ffd8b1', '#000075', '#808080', '#f0f0f0', '#000000')

# Reading the MLST Results
MLSTResults <- read.delim("MLSTResults.txt") %>% as_tibble()
MLSTResults <- MLSTResults %>% select(Sample, ST) %>% mutate(ST = gsub("\\*","", ST))
colnames(MLSTResults)[1] <- "Genome"

ann_colors <- list(ST = c("5" = colour[10], "7" = colour[15], "8" = colour[3], "9" = colour[9], "10" = colour[6], "11" = colour[17], "12" = colour[8], "39" = colour[13], "40" = colour[2], "41" = colour[16], "42" = colour[7], "43" = colour[4], "71" = colour[5], "88" = colour[12], "102" = colour[19], "Ancient" = colour[1], "NF" = colour[20], "Reference" = colour[22]))
#ann_colors = list(Pathovar = c("Ancient" = colour[1], "Commensal" = "#ffffff","EAEC" = colour[2], "EIEC" = "#22662A", "ETEC" = colour[16], "ExPEC" = colour[4],
#			       "Hybrid ExPEC/InPEC" = colour[3], "STEC" = colour[9],"Unknown" = colour[20]),
#		  #Host = c("Environment" = colour[10], "Food" = colour[11], "Human" = colour[12], "Livestock" = colour[13], "Porc" = colour[14], "Wild animal" = colour[15],"?" = colour[3]),
#		  ST = c("1429" = colour[6], "325" = colour[7], "3630" = colour[8], "399" = colour[12], "4995" = colour[18], "Other" = colour[20], "Ancient" = colour[1]),
#		  Ancient = c("Yes" = colour[22], "No" = colour[21]))
#
########### Pulling out the Phylogroups #######
#phylogroup <- read.csv("VirulenceHeatmapOlivier.csv") %>% as_tibble()
 
# What if we want the sources instead?
#olivierTable <- as_tibble(read.delim("PhyloMetadata.tab", header = T))
#olivierTable$Assembly.barcode <- sapply(olivierTable$Assembly.barcode, function(x){ifelse(grepl("^ETEC",x), gsub("ETEC", "ETEC ",x),x)})
#colnames(olivierTable)[1] <- "Genome"
#
#olivierTable$Pathovar <- sapply(olivierTable$Source, function(x){ifelse(grepl("Blood", x), "ExPEC", ifelse(grepl("Faeces",x), "Commensal",
#													 ifelse(grepl("Mummy",x), "Ancient", "Unknown")))})
#olivierTable$Pathovar <- sapply(1:length(olivierTable$ST), function(x){ifelse(grepl("ETEC", olivierTable$ST[x]), "ETEC",olivierTable$Pathovar[x])})
#olivierTable$Source <- sapply(olivierTable$Source, function(x){ifelse(grepl("animal|pig|,|Swine|Bovine|scrofa", x), "Animal",
#									    ifelse(grepl("Homo|ND|Blood|Faeces",x), "Human",
#										   ifelse(grepl("Cur|Cor|Ground", x), "Food", 
#											  ifelse(grepl("^$",x),"Unknown",x))))})
#
#olivierTable$Genome  <- sapply(olivierTable$Genome, function(x){ifelse(grepl("^ESC_", x), paste0(x,"_AS"),x)})
#olivierTable$Source  <- sapply(olivierTable$Source, function(x){ifelse(grepl("Human ", x),"Human",x)})
#olivierTable$Source  <- sapply(olivierTable$Source, function(x){ifelse(grepl("Mummy ", x),"Ancient",x)})

########### If we're cutting the heatmap in half #######
#filteredPhylogroup <- phylogroup %>% filter(Genome %in% strainList, !grepl("Shigella|Ancient|O104", Genome))
#preparation <- filteredPhylogroup %>% count(Phylogroup) %>% full_join(filteredPhylogroup, by = "Phylogroup") %>% mutate(Weight = 102/n)
#strainListV2  <- c(sample(preparation$Genome, size = nrow(preparation) * 0.5, prob = preparation$Weight), "O104")
#
#ind <- match(strainList, phylogroup$Genome)
#phylogroupHeat <- phylogroup$Phylogroup[ind]
#phylogroupHeat[length(phylogroupHeat)] <- "Ancient"
#
########### Pulling out the Pathovars #######
#pathovars <- as_tibble(read.csv("~/Downloads/PATRIC_genome.csv"))[,c("Genome.ID", "Pathovar")]
#colnames(pathovars)[1] <- "Genome"

############ Let's get the Depths ###########
op <- pboptions(type = "timer")
ncores = 8
depthDf <- do.call(bind_rows, pblapply("JessSampleTrimmed.tab.gz", cl = ncores,function(f){
		tmp <- as_tibble(read.delim(f, header = F, col.names = c("Gene", "Pos", "Coverage")))
		tmp$Genome <- gsub(".*/", "", gsub("\\..*|_genomic","",f)) 	
		tmp <- tmp %>% group_by(Genome, Gene) %>%
			summarize(MeanCoverage = mean(Coverage), sdCoverage = sd(Coverage), PercentCoverage = sum(Coverage > 0)/length(Coverage), .groups = "drop") %>%
			mutate(CV = sdCoverage/MeanCoverage)
		return(tmp)
	}))

# Getting the GC Content of the Genes
gcContent <- read.table("TrimmedPanGenomeGC.tab", header = F, col.names = c("Gene", "GC"), sep = "\t") %>% as_tibble() %>% mutate(GC = GC/100)
depthDf <- depthDf %>% left_join(gcContent)

# Now to correct for GC Bias
GCBiased <- depthDf %>% filter(MeanCoverage > 0) %>% select(Gene,MeanCoverage, GC)
colnames(GCBiased) <- c("Gene","Mean", "GCContent")
GCCorrected <- CorrectingGC(GCBiased)
depthDf <- depthDf %>% full_join(GCCorrected[[2]] %>% select(Gene, GCCorrected)) %>% mutate(GCCorrected = replace(GCCorrected, is.na(GCCorrected), 0))
depthDf$Genome <- "JessSample"

############# Core Gene Presence #####
roaryOutput <- as_tibble(read.delim("TrimmmedPresenceAbsence.tab"))
colnames(roaryOutput)[2:ncol(roaryOutput)] <- gsub("X|\\.scaffold|\\.genome|\\.result|_genomic", "", colnames(roaryOutput)[2:ncol(roaryOutput)])
colnames(roaryOutput)[2:ncol(roaryOutput)] <- gsub("\\.", "-", colnames(roaryOutput)[2:ncol(roaryOutput)])
# Finding the Core Genome
coreGenes <- rowSums(roaryOutput[,-1]) >= floor(0.95 * ncol(roaryOutput))
#coreGenes <- rowSums(ecoliOnly[,-1]) >= floor(0.99 * ncol(ecoliOnly))
coreGenes <- roaryOutput$Gene[coreGenes]

########### Some basic Coverage Data  ###############
ancientData <- depthDf %>% filter(GCCorrected >= 10, CV <=1) %>%
	pivot_wider(names_from = c(Genome), values_from=GCCorrected) %>% select(-c(MeanCoverage,GC,sdCoverage, PercentCoverage, CV)) #%>%
	#mutate(Genome = "10x - CV Filter")

ancientData %>% mutate(Status = ifelse(Gene %in% coreGenes, "Core", "Accessory")) %>% ggplot(aes(x = JessSample, fill = Status)) + 
	geom_histogram(position = "identity", alpha = 0.75, colour = "black", binwidth = 1) +
	scale_fill_manual(values = c(Accessory = "#007dba",Core = "#f8333c")) + theme_bw() +
	xlab("Mean Read Depth") + ylab("Genes") + theme(legend.position = "bottom") + 
	scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) + scale_y_continuous(breaks = scales::pretty_breaks(n = 10), limits = c(0, 750))
ggsave("GenePresenceHistogramTrimmed.pdf", width = 6, height = 4)

#####################################################################
# Gene Presence table
AllPA <- roaryOutput %>% left_join(ancientData) %>% mutate(JessSample = ifelse(is.na(JessSample)|JessSample == 0, 0,1))

# Now to make a boxplot to compare our gene to everyone else
tmp <- AllPA %>% select(-Gene) %>% summarize_all(sum) %>% t() %>% as.data.frame()
tmp$Genome <- rownames(tmp)
geneCounts <- tmp %>% as_tibble() 
colnames(geneCounts)[1] <- "Genes"
geneCounts %>% filter(Genome != "JessSample") %>% ggplot(aes(y = Genes)) + geom_boxplot() + theme_bw() +
	geom_point(data = geneCounts %>% filter(Genome == "JessSample"), aes(colour = Genome, x = 0)) +
	scale_colour_manual(values = colour) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), legend.position = "bottom")
ggsave("GeneCountsBmelTrimmed.pdf", width = 6, height = 9)

perCovMean <- depthDf %>%
       	filter(MeanCoverage > 0) %>% ggplot(aes(x= PercentCoverage,y = MeanCoverage)) + geom_point() + theme_bw() + scale_x_continuous(breaks = pretty_breaks(n = 10)) +
	scale_y_continuous(breaks = pretty_breaks(n = 10)) + ggtitle("Percent Coverage vs Mean Gene Coverage")
perCovMean <- ggMarginal(perCovMean, type = "densigram")
CVMean <- depthDf %>% filter(MeanCoverage > 0) %>% ggplot(aes(x= CV,y = MeanCoverage)) + geom_point() + theme_bw()+ scale_x_continuous(breaks = pretty_breaks(n = 10)) +
       	scale_y_continuous(breaks = pretty_breaks(n = 10))+ ggtitle("Coefficient of Variation vs Mean Gene Coverage")
CVMean <- ggMarginal(CVMean, type = "densigram")

pdf("~/CoverageMetricsTrimmed.pdf", width = 9, height = 9)
perCovMean
plot.new()
CVMean
dev.off()
#####################################################################

# Going to do this in a more sane way
coreDf <- roaryOutput %>% filter(Gene %in% coreGenes)
tmp <- colnames(coreDf)
tmp <- sapply(tmp, function(x){
	       ifelse(x %in% etecTranslations$GenomeID, etecList[x], x)
		  }) %>% unlist() %>% gsub(pattern = "\\.result.*|_genomic|\\.scaffold", replacement = "") %>%
	gsub(pattern = "\\.",replacement = "-")
colnames(coreDf) <- tmp

ancientData <- depthDf %>% filter(Genome %in%  c("JessSample"), GCCorrected >= 10, PercentCoverage >= 0.9) %>% mutate(Status = ifelse(Gene %in% coreGenes, "Core", "Accessory")) %>% filter(Status == "Core") %>%
	pivot_wider(names_from = c(Genome), values_from=GCCorrected) %>% select(-c(MeanCoverage, GC,sdCoverage, PercentCoverage, CV, Status))

#KosticData <- depthDf %>% filter(Genome %in% c("Lib4_3_bin","Lib4_7_bin","Lib4_8_bin"), MeanCoverage >= 1, PercentCoverage >= 0.9) %>%
#       	mutate(Status = ifelse(Gene %in% coreGenes, "Core", "Accessory")) %>%
#       	filter(Status == "Core") %>%
#	pivot_wider(names_from = c(Genome), values_from=MeanCoverage) %>%
#       	select(-c(sdCoverage, PercentCoverage, CV, Status))
#
#KosticData <- KosticData %>% group_by(Gene) %>% summarize_all(list(~ sum(., na.rm = T))) %>% arrange(Gene)
#
#KosticData %>% select(-Gene) %>% summarize_all(list(~ ifelse(. > 0, 1,.))) %>% summarize_all(sum)
#
#ancientData <- ancientData %>% full_join(KosticData) %>%
#	group_by(Gene) %>% summarize_all(list(~ sum(., na.rm = T))) %>%
#	ungroup() %>%
#       	mutate(KaeroEcoli = ifelse(is.na(KaeroEcoli) | KaeroEcoli == 0, 0, 1),
#		Lib4_3_bin = ifelse(is.na(Lib4_3_bin)| Lib4_3_bin == 0, 0, 1),
#		Lib4_7_bin = ifelse(is.na(Lib4_7_bin)| Lib4_7_bin == 0, 0, 1),
#		Lib4_8_bin = ifelse(is.na(Lib4_8_bin)| Lib4_8_bin == 0, 0, 1))# %>%
##	select(-Gene) %>% summarize_all(list(~ ifelse(. > 0, 1,.))) %>% summarize_all(sum)

coreData <- coreDf %>% left_join(ancientData) %>%
       	mutate(JessSample = ifelse(is.na(JessSample) | JessSample == 0, 0, 1))

# Here, we're trying to plot the core P/A data as a heatmap and/or a hierachical cluster
#coreData <- depthDf %>% filter(Genome %in% strainListOlivier, MeanCoverage >= 1, PercentCoverage >= 0.9) %>% mutate(Status = ifelse(Gene %in% coreGenes, "Core", "Accessory")) %>% filter(Status == "Core") %>%
#	pivot_wider(names_from = c(Genome), values_from=MeanCoverage) %>% select(-c(sdCoverage, PercentCoverage, CV, Status))

#ancientData <- depthDf %>% filter(Genome %in% strainListFragmented, MeanCoverage >= 10, PercentCoverage >= 0.9) %>% mutate(Status = ifelse(Gene %in% coreGenes, "Core", "Accessory")) %>% filter(Status == "Core") %>%
#	pivot_wider(names_from = c(Genome), values_from=MeanCoverage) %>% select(-c(sdCoverage, PercentCoverage, CV, Status))
#coreData <- coreData %>% full_join(ancientData)
#
#coreData[,-1] <- sapply(coreData[,-1], function(x){ifelse(is.na(x), 0, x)})
#coreData <- coreData %>% group_by(Gene) %>% summarize_all(sum) %>% group_by(Gene) %>% mutate_all(as.logical)

# Preparing the data
coreData <- as.data.frame(coreData)
rownames(coreData) <- coreData[,1]
coreData <- coreData[,-1]
genomes <- colnames(coreData)
coreData <- apply(coreData, MARGIN = 1,FUN = as.numeric)
rownames(coreData) <- genomes

coreDataSub <- coreData
maxCount <- coreDataSub %>% colSums() %>% max()
ind <- coreDataSub %>% colSums() == maxCount
coreDataSub <- coreDataSub[,!ind]

# Getting the heatmap sorted
#rownames(coreDataSub)[which(rownames(coreDataSub) %in% c("KaeroEcoli","Lib4_3_bin","Lib4_7_bin","Lib4_8_bin"))] <- c("AncientEcoli", "Zape2.1", "Zape2.2", "Zape3")
rownamesHeatmap <- ifelse(rownames(coreDataSub) == "AncientEcoli", "AncientEcoli","")
olivierTableHeatmap <- phylogroup[,c(4,2)] %>% as.data.frame()
olivierTableHeatmap$Pathovar <- sapply(olivierTableHeatmap$Pathovar, function(x){ifelse(x == "?", "Unknown",x)})
rownames(olivierTableHeatmap) <- phylogroup$Genome
olivierTableHeatmap <- olivierTableHeatmap[rownames(coreDataSub),]
#olivierTableHeatmap <- olivierTableHeatmap[complete.cases(olivierTableHeatmap),] %>%
#       	bind_rows(data.frame(row.names = c("Zape2.1","Zape2.2","Zape3"),ST = c("Kostic", "Kostic", "Kostic"), Pathovar = "Kostic"))
olivierTableHeatmap["AncientEcoli","ST"] <- "Ancient"

rownamesHeatmap <- ifelse(rownames(coreDataSub) == "AncientEcoli", "AncientEcoli","")
#rownames(olivierTableHeatmap) <- phylogroup$Genome
#rownames(olivierTableHeatmap) <- olivierTable$Genome

# Dendrogram and order
specDist <- dist(coreDataSub, method = "binary")
clusteredSpec <- hclust(specDist, method = "ward.D2") # 
rownameOrder <- clusteredSpec$order

treeColours <- as.list(ann_colors[[2]])[olivierTableHeatmap[rownameOrder,]$ST] %>% unlist()
treeColours <- gsub("#ffffff", "#000000", treeColours)
dend <- as.dendrogram(clusteredSpec)
labels_colors(dend) <- treeColours

pdf("~/CVFigures/TMPCoreTree.pdf", width = 24, height = 9)
#pdf("~/Documents/University/EcoliPaperV2/AdditionalFiles/PartsofFigures/SpeciesCoreTree.pdf", width = 6, height = 8)
plot(dend, horiz = F)
legend("topright", legend = names(ann_colors[[2]]), fill = unlist(ann_colors[[2]]), title = "ST")
dev.off()

# Now to cluster the trees
geneDist <- dist(coreDataSub)
geneClust <- hclust(geneDist, method = "ward.D2")

sil <- sapply(2:15, function(x){
	       tmp <- cutree(geneClust, k = x)
   	       summary(silhouette(tmp, geneDist))$avg.width
	 })

plot(y = sil, x = 2:15, type = "b")

small <- cutree(geneClust, k = 5)
clusterBorder <- colour[small] %>% unique()

pdf("~/KosticHierCoreClust.pdf", width = 12, height = 9)
plot(geneClust, hang = -1)
#rect.hclust(geneClust, k = 4, border = clusterBorder)
dev.off()

pheatmap(coreDataSub[rownameOrder,], annotation_row = olivierTableHeatmap,
	 color = c("#007dba", "#f8333c"), clustering_method = "ward.D2", legend = F, show_colnames = F,
	 gaps_row = c(which(rownames(coreDataSub)[rownameOrder] == "AncientEcoli"),which(rownames(coreDataSub)[rownameOrder] == "AncientEcoli") - 1),cluster_rows = F,
	 labels_row = rownamesHeatmap[rownameOrder], fontsize_row = 5, 
	 width = 12, height = 8, border_color = NA, annotation_colors = ann_colors, filename = "~/CVFigures/CorePACV.png")#, main = "Core Presence/Absence")
	 #filename = "~/Documents/University/EcoliPaperV2/AdditionalFiles/PartsofFigures/TestCore.png", width = 12, height = 8, border_color = NA, annotation_colors = ann_colors)#, main = "Core Presence/Absence")
#pheatmap(coreDataSub[rownameOrder,], annotation_row = olivierTableHeatmap,
#	 color = c("#007dba", "#f8333c"), clustering_method = "ward.D2", legend = F, show_colnames = F,
#	 gaps_row = c(which(rownames(coreDataSub)[rownameOrder] == "AncientEcoli"),which(rownames(coreDataSub)[rownameOrder] == "AncientEcoli") - 1),cluster_rows = F,
#	 cutree_cols = 4, labels_row = rownamesHeatmap[rownameOrder], fontsize_row = 5, 
#	 filename = "~/TestCore.png", width = 12, height = 8, border_color = NA, annotation_colors = ann_colors)#, main = "Core Presence/Absence")

# Core data counts for phylogeny
counts <- coreData %>% rowSums()
counts <- tibble(Genome = names(counts), Count = counts)
write.table(counts, sep = "\t", "../CoreGeneCounts.tab")

########## Accessory Genome Analysis ########
accessDf <- roaryOutput %>% filter(!(Gene %in% coreGenes))
#colnames(accessDf) %>% filter(!(Genome %in% c("IAI1", "ESC_NB8751AA_AS", "ATCC11231", "ESC_LB2165AA_AS")))
ancientData <- depthDf %>% filter(GCCorrected >= 10, CV <= 1) %>% mutate(Status = ifelse(Gene %in% coreGenes, "Core", "Accessory")) %>% filter(Status == "Accessory") %>%
	pivot_wider(names_from = c(Genome), values_from=GCCorrected) %>% select(-c(MeanCoverage,GC,sdCoverage, PercentCoverage, CV, Status))

ancientData <- ancientData %>%
	group_by(Gene) %>% summarize_all(list(~ sum(., na.rm = T))) %>%
	ungroup() %>%
       	mutate(JessSample = ifelse(is.na(JessSample) | JessSample == 0, 0, 1))# %>%

accessData <- accessDf %>% left_join(ancientData) %>%# mutate(JessSample = ifelse(is.na(JessSample), 0, 1))
       	mutate(JessSample = ifelse(is.na(JessSample) | JessSample == 0, 0, 1))

# Now, we'll be doing the same thing for the accessory genome.  Firstly, however, is going to be the removal genes which are barely present.  We'll do a separate one with unique genes
#accessData <- depthDf %>% filter(Genome %in% strainListOlivier, MeanCoverage >= 1) %>% mutate(Status = ifelse(Gene %in% coreGenes, "Core", "Accessory")) %>% filter(Status == "Accessory") %>%
	#pivot_wider(names_from = c(Genome), values_from=MeanCoverage) %>% select(-c(sdCoverage, PercentCoverage, CV, Status))
#ancientData <- depthDf %>% filter(Genome %in% strainListFragmented, MeanCoverage >= 10) %>% mutate(Status = ifelse(Gene %in% coreGenes, "Core", "Accessory")) %>% filter(Status == "Accessory") %>%
#	pivot_wider(names_from = c(Genome), values_from=MeanCoverage) %>% select(-c(sdCoverage, PercentCoverage, CV, Status))
#
#accessData <- accessData %>% full_join(ancientData)

accessData[,-1] <- sapply(accessData[,-1], function(x){ifelse(is.na(x), 0, x)})
accessData <- accessData %>% group_by(Gene) %>% summarize_all(sum) %>% group_by(Gene) %>% mutate_all(as.logical)

# Preparing the data
accessData <- as.data.frame(accessData)
rownames(accessData) <- accessData[,1]
accessData <- accessData[,-1]
genomes <- colnames(accessData)
accessData <- apply(accessData, MARGIN = 1,FUN = as.numeric)
rownames(accessData) <- genomes

# We'll be removing the unique genes now
accessDataSub <- accessData
maxCount <- accessDataSub %>% colSums() %>% max()
ind <- accessDataSub %>% colSums() < 5
accessDataSub <- accessDataSub[,!ind]

# Getting the Metadata Setup
#rownames(accessDataSub)[which(rownames(accessDataSub) %in% "KaeroEcoli")] <-"AncientEcoli"
#rownamesHeatmap <- ifelse(rownames(accessDataSub) == "AncientEcoli", "AncientEcoli","")
#olivierTableHeatmap <- phylogroup[,c(4,2)] %>% as.data.frame()
#olivierTableHeatmap$Pathovar <- sapply(olivierTableHeatmap$Pathovar, function(x){ifelse(x == "?", "Unknown",x)})
#rownames(olivierTableHeatmap) <- phylogroup$Genome
#olivierTableHeatmap <- olivierTableHeatmap[rownames(accessDataSub),]
#olivierTableHeatmap["AncientEcoli","ST"] <- "Ancient"

# Dendrogram and order
specDist <- dist(accessDataSub, method = "binary")
clusteredSpec <- hclust(specDist, method = "ward.D2") # 
rownameOrder <- clusteredSpec$order

treeColours <- ifelse(rownames(accessDataSub)[rownameOrder] == "JessSample", "red", "black")
dend <- as.dendrogram(clusteredSpec)
labels_colors(dend) <- treeColours

pdf("AccessTreeTrimmed.pdf", width = 48, height = 18)
#pdf("~/Documents/University/EcoliPaperV2/AdditionalFiles/PartsofFigures/SpeciesAccessTree.pdf", width = 6, height = 8)
plot(dend, horiz = F, main = "Accessory Genome Clustering")
dev.off()

# Now to cluster the trees
#geneDist <- dist(t(accessDataSub))
geneDist <- dist(accessDataSub, method = "euclidean")
geneClust <- hclust(geneDist, method = "ward.D2")
######################################
# If I want to use a MDS plot
fit <- cmdscale(geneDist, eig = T, k = 4, add =T)
coord <- fit$points %>% as_tibble()
coord$Genome <- rownames(fit$points)
contrib <- fit$eig/sum(fit$eig) * 100
coord <- coord %>% left_join(MLSTResults)
coord$ST[which(coord$Genome == "JessSample")] <- "Ancient" 

# Quickly filtering the list so that only the STs present are used
ann_colors$ST <- ann_colors$ST[unique(coord$ST)]
p1 <- coord %>%
       	ggplot(aes(x = V1, y = V2, label = Genome, colour = ST)) +
	geom_hline(yintercept=0, lty = 2, colour = "grey90") + 
	geom_vline(xintercept=0, lty = 2, colour = "grey90") +
	geom_point() +
	xlab(bquote("PCoA 1 ("~.(round(contrib[1],2))~"%)")) +
	ylab(bquote("PCoA 2 ("~.(round(contrib[2],2))~"%)")) +
	scale_colour_manual(values = ann_colors$ST) +
	#geom_text_repel(show.legend = F) +
	theme_classic()
p1

p2 <- coord %>%
       	ggplot(aes(x = V3, y = V4, label = Genome, colour = ST)) +
	geom_hline(yintercept=0, lty = 2, colour = "grey90") + 
	geom_vline(xintercept=0, lty = 2, colour = "grey90") +
	geom_point() +
	xlab(bquote("PCoA 3 ("~.(round(contrib[3],2))~"%)")) +
	ylab(bquote("PCoA 4 ("~.(round(contrib[4],2))~"%)")) +
	scale_colour_manual(values = ann_colors$ST) +
	#geom_text_repel(show.legend = F) +
	theme_classic()
p2

ggarrange(p1,p2, legend = "bottom", align = "hv", ncol = 1, common.legend = T)
ggsave("PCoA_AccessoryTrimmed.pdf", width = 9, height = 6)

# Clustering based on PCoA Coordinates
sil <- sapply(2:6, function(i){clara(coord[,-c(5,6)], i)$silinfo$avg.width})
plot(2:6, sil, type = "b")

clustered <- clara(coord[,-c(5,6,6)], 3)

coord$clusters <- factor(clustered$clustering)

p1 <- coord %>%
       	ggplot(aes(x = V1, y = V2, label = Genome, colour = ST, shape = clusters, group = clusters)) +
	geom_hline(yintercept=0, lty = 2, colour = "grey90") + 
	geom_vline(xintercept=0, lty = 2, colour = "grey90") +
	geom_point() +
	stat_ellipse() +
	xlab(bquote("PCoA 1 ("~.(round(contrib[1],2))~"%)")) +
	ylab(bquote("PCoA 2 ("~.(round(contrib[2],2))~"%)")) +
	scale_colour_manual(values = ann_colors$ST) +
	#geom_text_repel(show.legend = F) +
	theme_classic()
p1

p2 <- coord %>%
       	ggplot(aes(x = V3, y = V4, label = Genome, colour = ST, shape = clusters, group = clusters)) +
	geom_hline(yintercept=0, lty = 2, colour = "grey90") + 
	geom_vline(xintercept=0, lty = 2, colour = "grey90") +
	geom_point() +
	stat_ellipse() +
	xlab(bquote("PCoA 3 ("~.(round(contrib[3],2))~"%)")) +
	ylab(bquote("PCoA 4 ("~.(round(contrib[4],2))~"%)")) +
	scale_colour_manual(values = ann_colors$ST) +
	#geom_text_repel(show.legend = F) +
	theme_classic()
p2

ggarrange(p1,p2, legend = "bottom", align = "hv", ncol = 1, common.legend = T)
ggsave("PCoA_Accessory_Cluster.pdf", width = 9, height = 6)
ggsave(p1,file = "PCoA_Accessory_ClusterAxes12.pdf", width = 9, height = 6)

###########################################################################
# Now, I'll redo the analysis, however, only using genomes from cluster 1 #
###########################################################################
Clust1Genomes <- coord %>% filter(clusters == 1) %>% pull(Genome)

accessDataClust1 <- accessData[Clust1Genomes,]
maxCount <- accessDataClust1 %>% colSums() %>% max()
ind <- accessDataClust1 %>% colSums() < 5
accessDataClust1 <- accessDataClust1[,!ind]

geneDist <- dist(accessDataClust1, method = "euclidean")
geneClust <- hclust(geneDist, method = "ward.D2")

# Making the MDS Plot
fit <- cmdscale(geneDist, eig = T, k = 4, add =T)
coord <- fit$points %>% as_tibble()
coord$Genome <- rownames(fit$points)
contrib <- fit$eig/sum(fit$eig) * 100
coord <- coord %>% left_join(MLSTResults)
coord$ST[which(coord$Genome == "JessSample")] <- "Ancient" 

p1 <- coord %>%
       	ggplot(aes(x = V1, y = V2, label = Genome, colour = ST)) +
	geom_hline(yintercept=0, lty = 2, colour = "grey90") + 
	geom_vline(xintercept=0, lty = 2, colour = "grey90") +
	geom_point() +
	xlab(bquote("PCoA 1 ("~.(round(contrib[1],2))~"%)")) +
	ylab(bquote("PCoA 2 ("~.(round(contrib[2],2))~"%)")) +
	scale_colour_manual(values = ann_colors) +
	#geom_text_repel(show.legend = F) +
	theme_classic()
p1

p2 <- coord %>%
       	ggplot(aes(x = V3, y = V4, label = Genome, colour = ST)) +
	geom_hline(yintercept=0, lty = 2, colour = "grey90") + 
	geom_vline(xintercept=0, lty = 2, colour = "grey90") +
	geom_point() +
	xlab(bquote("PCoA 3 ("~.(round(contrib[3],2))~"%)")) +
	ylab(bquote("PCoA 4 ("~.(round(contrib[4],2))~"%)")) +
	scale_colour_manual(values = ann_colors) +
	#geom_text_repel(show.legend = F) +
	theme_classic()
p2

ggarrange(p1,p2, legend = "bottom", align = "hv", ncol = 1, common.legend = T)
ggsave("AccessoryPCoA_Cluster1ONLY.pdf", width = 9, height = 6)

# Clustering based on PCoA Coordinates
sil <- sapply(2:15, function(i){clara(coord[,-c(5,6)], i)$silinfo$avg.width})
plot(2:15, sil, type = "b")

clustered <- clara(coord[,-c(5,6)], 3)
coord$clusters <- factor(clustered$clustering)

p1 <- coord %>%
       	ggplot(aes(x = V1, y = V2, label = Genome, colour = ST, shape = clusters, group = clusters)) +
	geom_hline(yintercept=0, lty = 2, colour = "grey90") + 
	geom_vline(xintercept=0, lty = 2, colour = "grey90") +
	geom_point() +
	stat_ellipse() +
	xlab(bquote("PCoA 1 ("~.(round(contrib[1],2))~"%)")) +
	ylab(bquote("PCoA 2 ("~.(round(contrib[2],2))~"%)")) +
	scale_colour_manual(values = ann_colors) +
	#geom_text_repel(show.legend = F) +
	theme_classic()
p1

p2 <- coord %>%
       	ggplot(aes(x = V3, y = V4, label = Genome, colour = ST, shape = clusters, group = clusters)) +
	geom_hline(yintercept=0, lty = 2, colour = "grey90") + 
	geom_vline(xintercept=0, lty = 2, colour = "grey90") +
	geom_point() +
	stat_ellipse() +
	xlab(bquote("PCoA 3 ("~.(round(contrib[3],2))~"%)")) +
	ylab(bquote("PCoA 4 ("~.(round(contrib[4],2))~"%)")) +
	scale_colour_manual(values = ann_colors) +
	#geom_text_repel(show.legend = F) +
	theme_classic()
p2

ggarrange(p1,p2, legend = "bottom", align = "hv", ncol = 1, common.legend = T)
ggsave("PCoA_Accessory_Cluster1SUBCLUSTERS3.pdf", width = 9, height = 6)
ggsave(p1,file = "PCoA_Accessory_ClusterAxes12.pdf", width = 9, height = 6)
#####################################
# Seeing what happens if I use a MCA
tmp <- data.frame(accessDataSub)
tmp <- sapply(tmp, as.logical)
MCAResults <- MCA(tmp, graph = F)
tmpHabillage <- ifelse(rownames(tmp) == "JessSample", "Ancient", "Modern") 
fviz_mca_ind(MCAResults, habillage = factor(tmpHabillage), geom = "point", palette= c("black", "red"), axes = c(1,2))

#####################################

sil <- sapply(2:15, function(x){
	       tmp <- cutree(geneClust, k = x)
   	       summary(silhouette(tmp, geneDist))$avg.width
	 })

plot(y = sil, x = 2:15, type = "b")

small <- cutree(geneClust, k = 4)
clusterBorder <- colour[small] %>% unique()

annotationData <- MLSTResults %>% bind_rows(data.frame(Genome = "JessSamples",ST = "Ancient")) %>% as.data.frame()
annotationData2 <- as.data.frame(annotationData[,-1])
rownames(annotationData2) <- annotationData$Genome
annotationData <- annotationData2
rm(annotationData2)
colnames(annotationData) <- "ST"
#annotationData$clusters <- as.factor(annotationData$clusters)
#ann_colors <- list(Colour = c("Modern" = "#ffffff", "Ancient" = colour[1]), clusters = c("1" = colour[6], "2" = colour[7], "3" = colour[8], "4" = colour[9]))
ann_colors_heatmap <- list(ST = unlist(ann_colors))
#ann_colors_heatmap <- list("ST" = ann_colors)

pheatmap(accessDataSub, annotation_row = annotationData,
	 color = c("#007dba", "#f8333c"), clustering_method = "ward.D2", legend = F, show_colnames = F,
	show_rownames = F, 
	 #border_color = NA, annotation_colors = ann_colors)#, main = "Accessory Presence/Absence")
	 filename = "AccessoryHeatmapTrimmed.pdf", width = 12, height = 8, border_color = NA, annotation_colors = ann_colors)#, main = "Accessory Presence/Absence")
