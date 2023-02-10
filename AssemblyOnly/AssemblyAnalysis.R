library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(ape)
library(scales)
library(ggExtra)
library(ggnewscale)
#library(ggbreak)
library(cluster)
library(parallel)
library(pbapply)
library(purrr)
library(pheatmap)
library(dendextend)
library(factoextra)
library(gtools)
library(xtable)
library(eulerr)
library(rjson)
## Functions
BlastParsing <- function(blastFile, ncores){
	tmp <- read.delim(blastFile, header = F,
			  col.names = c("Query", "Match", "PIdent", "Length", "Mismatches", "Gapopen", "QStart", "QEnd", "SStart", "SSend", "Evalue", "Bitscore", "Taxa")) %>% as_tibble()
	
	tmpList <- split(tmp, as.factor((tmp$Query)))
	
	op <- pboptions(type = "timer")
	FilteredAmbiguousGenes<- pblapply(cl = ncores,tmpList, function(x){ 
		       	filtResults <- x %>% filter(Evalue == min(Evalue))
	
			if(nrow(filtResults) > 1){
				filtResults <- filtResults %>% filter(Bitscore == max(Bitscore))
				if(nrow(filtResults) > 1){
					filtResults <- filtResults %>% filter(PIdent == max(PIdent))
					}
					if(nrow(filtResults) > 1){
						filtResults <- filtResults %>% filter(Mismatches == min(Mismatches))
						}
				}
			return(filtResults[1,])
	 }) %>% bind_rows() %>% filter(PIdent >= 90)
}

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
theme_set(theme_classic())

colour <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6',
'#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3',
'#808000', '#ffd8b1', '#000075', '#808080', '#f0f0f0', '#000000')

clusterColours <- c("Western Mediterranean" = "#4D9DE0", "Fertile Crescent" = "#E15554", "Africa/America" = "#E1BC29", "Indo-Pacific" = "#3BB273", "Russia" = "#7768AE")
tmp <- read.delim("CountriesGrouped.list", header = F, col.names = c("Country", "Colour"))
countryColours <- as.list(tmp$Colour)
names(countryColours) <- tmp$Country
rm(tmp)

# Reading the MLST Results
MLSTResults <- read.delim("MLSTResults.txt") %>% as_tibble()
MLSTResults <- MLSTResults %>% select(Sample, ST) %>% mutate(ST = gsub("\\*| ","", ST))
colnames(MLSTResults)[1] <- "Genome"

# Reading the metadata table
metaData <- read.delim("MetadataAll.tab", header =T) %>% as_tibble
metaData <- metaData[-266,] # Have duplicate....

ann_colors <- list(ST = c("5" = colour[3], "7" = colour[15], "8" = colour[10], "9" = colour[9], "10" = colour[6], "11" = colour[17], "12" = colour[8], "39" = colour[13], "40" = colour[2], "41" = colour[16], "42" = colour[7], "43" = colour[4], "71" = colour[11], "88" = colour[12], "102" = colour[19],"Ancient" = colour[1], "NF" = colour[20], "NIPH" = colour[22]),
	Norway = c("Normal" = "#00205B" , "Odd" = "#BA0C2F", "Other" = "#FFFFFF"))

############# Mean Assembly Coverage ###########
assemblyMapping <- read.delim(file = "SummarizedDepths.tab", header = F, col.names = c("Source", "CHROM", "Mean", "SD", "SE", "CV", "PCov")) %>% as_tibble()

gcContent <- read.delim("BrucellaContigsFilteredStats.tab", header = F, col.names = c("CHROM", "Length", "GCContent"))

assemblyMappingCorrected <- assemblyMapping %>% left_join(gcContent) %>% CorrectingGC()
assemblyMapping <- assemblyMappingCorrected$UpdatedDepths

summarizedResults <- assemblyMapping %>% 
	summarize(MeanCov = mean(GCCorrected), SD = sd(GCCorrected), SE = qnorm(0.975) * SD/sqrt(nrow(assemblyMapping)), Lo = MeanCov - SE, Hi = MeanCov + SE) %>%
       	as.data.frame()

assemblyMapping %>% mutate(Length = as.numeric(gsub(".*length_|_.*","", CHROM))) %>% arrange(-Mean) 

histPlot <- assemblyMapping %>% ggplot(aes(x = GCCorrected)) + 
	geom_rect(inherit.aes = F,data = summarizedResults,
		  aes(ymin = -Inf, ymax = Inf, xmin = MeanCov - 2 * SD, xmax =  MeanCov + 2 * SD), fill = "black", alpha = 0.5) +
	geom_histogram(fill = colour[1], colour = "black") +
	geom_vline(xintercept = summarizedResults$MeanCov, lty = 2) +
	xlab("Mean Read Depth") + ylab("Contigs")
	#scale_x_break(breaks = c(35,45), space = 1)

GCScatter <- assemblyMapping %>% select(CHROM, Mean, GCCorrected, GCContent) %>% pivot_longer(-c(CHROM,GCContent), names_to = "GCCorrection", values_to = "Coverage") %>% 
	mutate(GCCorrection = ifelse(GCCorrection == "Mean", "Before", "After"))

scatterPlot <- GCScatter %>% ggplot(aes(x = GCContent/100, y = Coverage, colour = GCCorrection)) +
	geom_point() +
	geom_smooth(method = "lm") +
	scale_y_continuous(breaks = breaks_pretty(10)) +
	scale_colour_manual(values = c("Before" = colour[4], "After" = colour[1]), "GC Correction") +
	ylab("Mean Read Coverage") +
	xlab("GC Content")

bot <- ggarrange(scatterPlot, plot.new(), nrow = 1, common.legend = T, legend = "bottom", labels = c("b","c"))
ggarrange(histPlot, bot, ncol = 1, labels = c("A",""))
ggsave("Figure2.pdf", width = 9, height = 6)

############# Core Gene Presence #####
roaryOutput <- as_tibble(read.delim("gene_presence_absence.Rtab"))
colnames(roaryOutput)[2:ncol(roaryOutput)] <- gsub("^X|\\.scaffold|\\.genome|\\.result|_genomic", "", colnames(roaryOutput)[2:ncol(roaryOutput)])
colnames(roaryOutput)[2:ncol(roaryOutput)] <- gsub("\\.", "-", colnames(roaryOutput)[2:ncol(roaryOutput)])
# Finding the Core Genome
coreGenes <- rowSums(roaryOutput[,-1]) >= floor(0.99 * ncol(roaryOutput))
coreGenes <- roaryOutput$Gene[coreGenes]

# Getting as much information as possible on the missing core genes
PABrancorsini <- roaryOutput[,c(1,169)]

coreGenesMissingBrancorsini <- PABrancorsini %>% left_join(read.delim("../BlastResults/Pangenome.tab"), by = c("Gene" = "Query")) %>% rename(COGGene = Gene.y) %>% 
	filter(Brancorsini == 0, Gene %in% coreGenes)
coreGenesMissingBrancorsini$Gene <- apply(coreGenesMissingBrancorsini, MARGIN = 1, function(x){
	      if(grepl("group_",x[1]) & !is.na(x[8]) & nchar(x[8]) > 0){
		      return(x[8])
	      }else{
		      return(x[1])
	      }
	})

coreGenesMissingBrancorsini$Gene %>% write.table("CoreMissing.list", col.names = F, row.names =F, quote =F)

rm(coreGenesMissingGeridu, coreGenesMissingBrancorsini)

# Getting the Gene Count Boxplot
tmp <- roaryOutput %>% select(-Gene) %>% summarize_all(sum) %>% t() %>% as.data.frame()
tmp$Genome <- rownames(tmp)
geneCounts <- tmp %>% as_tibble() 
geneCounts <- geneCounts %>% left_join(MLSTResults)
colnames(geneCounts)[1] <- "Genes"
geneCounts %>% filter(grepl("Brancorsini|Geridu", Genome))
geneCounts %>% filter(!grepl("Brancorsini|Geridu", Genome)) %>% summarize(mean(Genes))
geneCounts %>% filter(!grepl("Brancorsini|Geridu", Genome)) %>% ggplot(aes(y = Genes, x ="", color = ST)) +
       	geom_boxplot(colour = "black") + 
	#geom_jitter(height = 0) +
	geom_point(data = geneCounts %>% filter(grepl("Brancorsini|Geridu", Genome)), aes(colour = Genome, x = "")) +
	scale_colour_manual(values = ann_colors$ST[names(ann_colors$ST) %in% c("Brancorsini", "Geridu")]) + 
	theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), legend.position = "right") +
	guides(colour = guide_legend(title = "Sample"))
ggsave("GeneCountsAssembly.pdf", width = 6, height = 9)

# Any genes unique to Brancorsini?
roaryDF <- roaryOutput %>% as.data.frame()
rownames(roaryDF) <- roaryDF$Gene
roaryDF <- roaryDF[,-1]
genomeCount <- rowSums(roaryDF)
roaryOutput %>% select(Gene, Brancorsini) %>% filter(Gene %in% names(genomeCount)[genomeCount == 1], Brancorsini > 0)

##############################
### Core Gene Scatter Plot ###
##############################
coreDf <- roaryOutput %>% filter(Gene %in% coreGenes)

# Preparing the data
coreData <- as.data.frame(coreDf)
rownames(coreData) <- coreData[,1]
coreData <- coreData[,-1]
genomes <- colnames(coreData)
coreData <- apply(coreData, MARGIN = 1,FUN = as.numeric)
rownames(coreData) <- genomes

coreDataSub <- coreData
maxCount <- coreDataSub %>% colSums() %>% max()
ind <- coreDataSub %>% colSums() == maxCount
coreDataSub <- coreDataSub[,!ind]

######################################
# If I want to use a MDS plot
geneDist <- dist(coreDataSub, method = "binary")

fit <- cmdscale(geneDist, eig = T, k = 4, add =T)
coord <- fit$points %>% as_tibble()
coord$Genome <- rownames(fit$points)
contrib <- fit$eig/sum(fit$eig) * 100
coord <- coord %>% left_join(MLSTResults)
coord$ST[which(coord$Genome == "Brancorsini")] <- c("Ancient")

# Quickly filtering the list so that only the STs present are used
coreColours <- ann_colors$ST[names(ann_colors$ST) %in% coord$ST]
pCore <- coord %>%
       	ggplot(aes(x = V1, y = V2, label = Genome, colour = ST)) +
	geom_hline(yintercept=0, lty = 2, colour = "grey90") + 
	geom_vline(xintercept=0, lty = 2, colour = "grey90") +
	geom_point(alpha = 0.75) +
	xlab(bquote("PCoA 1 ("~.(round(contrib[1],2))~"%)")) +
	ylab(bquote("PCoA 2 ("~.(round(contrib[2],2))~"%)")) +
	scale_colour_manual(values = coreColours) +
	guides(colour = guide_legend(nrow = 2)) +
	#geom_text_repel(show.legend = F) +
	theme_classic() +
	theme(legend.position = "bottom") 
pCore
ggsave(pCore, file = "PCoACoreDec2022.pdf", width = 6, height = 4)

corePlotData <- coord %>% select(Genome, ST) %>% distinct()%>%
       	mutate(ST = replace(ST, grepl("^NIPH-*|NI_2007", Genome) & grepl("8",ST), "NIPH")) %>%
	mutate(ST = replace(ST, grepl("Brancorsini", Genome), "Ancient")) %>%
	as.data.frame()

annotationHeat <- as.data.frame(corePlotData$ST)
rownames(annotationHeat) <- corePlotData$Genome
colnames(annotationHeat) <- "ST"

ann_colorsHeat <- list(ST = ann_colors$ST)
pheatmap(coreDataSub, clustering_distance_rows = "binary",clustering_distance_cols = "binary", clustering_method = "ward.D2",
	 show_rownames = F, show_colnames = F, legend = F,
	 annotation_row = annotationHeat, annotation_colors = ann_colorsHeat, annotation_legend = T,
	 annotation_names_row = T, filename = "CoreHeatmap.pdf", width = 6, height = 6)

#################################
### Accessory Genome Analysis ###
#################################
accessDf <- roaryOutput %>% filter(!(Gene %in% coreGenes))

# Preparing the data
accessData <- as.data.frame(accessDf)
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

######################################
# If I want to use a MDS plot
geneDist <- dist(accessDataSub, method = "binary")

fit <- cmdscale(geneDist, eig = T, k = 4, add =T)
coord <- fit$points %>% as_tibble()
coord$Genome <- rownames(fit$points)
contrib <- fit$eig/sum(fit$eig) * 100
coord <- coord %>% left_join(MLSTResults)
coord$ST[coord$Genome %in% c("Brancorsini")] <- c("Ancient")

# Clustering based on PCoA Coordinates
fviz_nbclust(coord[,1:4],FUNcluster = clara, method = "silhouette")
fviz_nbclust(coord[,1:4],FUNcluster = clara, method = "wss")
#fviz_nbclust(coord[,1:4],FUNcluster = clara, method = "gap_stat")

clustered <- clara(coord[,-c(5,6)], 5)

coord$clusters <- clustered$clustering
coord  <- coord %>% mutate(clusters = ifelse(clusters == 1, "Western Mediterranean",
				 ifelse(clusters == 2, "Fertile Crescent",
				       	ifelse(clusters == 3, "Africa/America", ifelse(clusters == 4, "Indo-Pacific", "Russia")))))
#write.table(coord, file = "FullPhyloClustering.tab", sep = "\t", row.names = F, col.names = T)

accessPlotData <- coord %>% left_join(metaData %>% select(-ST), by = c("Genome" = "Sample")) %>%
       	#mutate(Norway3 = ifelse(grepl("^NIPH-*|NI_2007", Genome) & grepl("8",ST), "OddNorway", NA)) %>%
       	mutate(ST = replace(ST, grepl("^NIPH-*|NI_2007", Genome) & grepl("8",ST), "NIPH"))
coreColours <- ann_colors$ST[names(ann_colors$ST) %in% accessPlotData$ST]
coreColours <- coreColours[mixedsort(names(coreColours))]
accessPlotData$ST <- factor(accessPlotData$ST, levels = unique(mixedsort(accessPlotData$ST)))
accessPlotData <- accessPlotData %>% mutate(Type = ifelse(Genome == "Reference", "Reference", ifelse(Genome == "Brancorsini", "Ancient", "Modern"))) %>%
	mutate(Type = factor(Type, levels = c("Modern", "Ancient", "Reference")))

# Getting the labels for the ellipses
labelCoord <- accessPlotData %>% group_by(clusters) %>% summarize(V1 = mean(V1), V2 = mean(V2), V3 = mean(V3), V4 = mean(V4))

p1 <- accessPlotData %>%
       	ggplot(aes(x = V1, y = V2, label = Genome, colour = ST, shape = Type, group = clusters)) +#, group = clusters)) +
	geom_hline(yintercept=0, lty = 2, colour = "grey90") + 
	geom_vline(xintercept=0, lty = 2, colour = "grey90") +
	geom_point() +
	scale_colour_manual(values = coreColours, "Sequence Type") +
	guides(colour = guide_legend(nrow = 2, title.position = "top", title.hjust = 0.5)) +
	new_scale_colour() +
	stat_ellipse(mapping = aes(colour = clusters), show.legend = F) +
	scale_colour_manual("Accessory Gene Clustering", values = clusterColours) +
	#stat_ellipse(data = . %>% filter(!is.na(Norway3)),show.legend = F, colour = "black", lty = 2,mapping = aes(group = Norway3)) +
	xlab(bquote("PCoA 1 ("~.(round(contrib[1],2))~"%)")) +
	ylab(bquote("PCoA 2 ("~.(round(contrib[2],2))~"%)")) +
	geom_text_repel(data = labelCoord, aes(x = V1, y = V2, label = clusters, colour = clusters), inherit.aes = F, show.legend = F) +
	guides(colour = guide_legend(nrow = 2, title.position = "top", title.hjust = 0.5), shape = guide_legend(nrow = 1, title.position = "top", title.hjust = 0.5, title = "Sample")) +
	theme(legend.position = "bottom")
 #theme(legend.position = "bottom", axis.text.x = element_blank(), axis.title.x = element_blank())

ggsave(p1,file = "FixedAccessoryDecPCoA.pdf", width = 9, height = 6)

# Now to look at a heatmap
groupedCountries <- c("Afghanistan" = "Asia", "Argentina" = "Americas", "China" = "China", "Ethiopia" = "Africa", "Georgia" = "Asia", "India" = "India", "Iraq" = "Fertile Crescent","Israel" = "Fertile Crescent", "Italy" = "Italy", "Kuwait" = "Fertile Crescent", "Malta" = "Europe", "Nigeria" = "Africa", "Norway" = "Norway", "Protugal" = "Europe","Somalia" = "Africa", "Turkey" = "Fertile Crescent", "United Kingdom" = "Europe", "Zimbabwe" = "Africa", "Unknown" = "Unknown", "Albania" = "Europe","Bulgaria" = "Europe", "Canada" = "Americas", "Cyprus" = "Europe", "Iran" = "Fertile Crescent", "Jordan" = "Fertile Crescent", "Kosovo" = "Europe","Malaysia" = "Asia", "Pakistan" = "Asia", "Russia" = "Russia", "Saudi Arabia" = "Asia", "Sudan" = "Africa", "Syria" = "Fertile Crescent", "Thailand" = "Asia","Turkmenistan" = "Asia", "Egypt" = "Egypt", "France" = "Europe", "Morocco" = "Morocco", "United States" = "Americas", "Portugal" = "Europe")

annotationHeat <- accessPlotData %>% select(Genome, ST, clusters, Country) %>% distinct() %>%
	mutate(Country = replace(Country, is.na(Country), "Unknown"), Region = groupedCountries[Country]) %>%
	mutate(Region = factor(Region, levels = names(countryColours))) %>% select(-Country) %>%
       	as.data.frame()
#annotationHeat <- annotationHeat[-166,] # Duplicate being removed as it has no information

rownames(annotationHeat) <- annotationHeat$Genome
annotationHeat <- annotationHeat[,-1]
colnames(annotationHeat)[2] <- "PAV Clusters"

ann_colorsHeat <- list("PAV Clusters" = clusterColours, ST = ann_colors$ST, Region = unlist(countryColours))

pheatmap(accessDataSub, clustering_distance_rows = "binary",clustering_distance_cols = "binary", clustering_method = "ward.D2",
	 show_rownames = F, show_colnames = F, legend = F,
	 cutree_rows = 5,
	 annotation_row = annotationHeat, annotation_colors = ann_colorsHeat, annotation_legend = T,
	 annotation_names_row = T, filename = "AccessoryHeatmap.pdf", width = 12, height = 12)

# Let's only look at the third cluster now! Hopefully get names as well
#annotationHeat <- accessPlotData %>% select(Genome, ST, clusters) %>% filter(grepl("Afr", clusters)) %>% distinct() %>% as.data.frame()
#rownames(annotationHeat) <- annotationHeat$Genome
#annotationHeat <- annotationHeat[,-1]
#colnames(annotationHeat)[2] <- "Accessory Clusters"
#
#ann_colorsHeat <- list("Accessory Clusters" = clusterColours, ST = ann_colors$ST, Region = unlist(countryColours))
#africaOnly <- accessDataSub[rownames(annotationHeat),]
##africaDifferential <- africaOnly[,colSums(africaOnly) != 0]
#
#pheatmap(africaOnly, clustering_distance_rows = "binary",clustering_distance_cols = "binary", clustering_method = "ward.D2",
#	 show_rownames = F, show_colnames = F, legend = F,
#	 cutree_rows = 5,
#	 annotation_row = annotationHeat, annotation_colors = ann_colorsHeat, annotation_legend = T,
#	 annotation_names_row = F, filename = "AfricaOnly.pdf", width = 6, height = 6)
#
### What makes the three EMed clusters? ###
clusteredResults <- coord %>% mutate(clusters = as.character(clusters))

clusteredResults <- split(clusteredResults %>% pull(Genome), f = clusteredResults$clusters)

sharedGenesbyCluster <- lapply(clusteredResults, function(x){
	       tmp <- accessDf[,c("Gene",x)] %>% as.data.frame()
	       rownames(tmp) <- tmp$Gene
	       tmp <- tmp[,-1]
	       tmp <- rowSums(tmp)
	       return(names(tmp)[tmp > 1])
	 })
tmp <- sharedGenesbyCluster %>% unlist()

tmp <- tibble(Cluster = names(tmp),Gene = tmp) %>% mutate(Cluster = gsub("\\d*","", Cluster), TMP = T) %>%
       	pivot_wider(id_cols = Gene, names_from = Cluster, values_from = TMP, values_fill = F)
write.table(tmp, file = "AccessoryGeneLocation.tab", row.names = F, sep = "\t", quote = F)

# Making a quick table
geneLocation <- tmp
tmp <- geneLocation %>% filter(Gene %in% colnames(accessDataSub))

geneInfo <- tmp %>% mutate(ClustersPresent = rowSums(tmp[,-1])) %>%
		filter(ClustersPresent == 1 | ClustersPresent == 4)
       
geneInfo <- geneInfo %>% mutate(Unique = apply(geneInfo, 1, function(data){
					   if(data[7] == 1){
						   return(names(which(data == T)))
					   }else{
						   return(NA)
					   }
				 })) %>%
		mutate(Absent = apply(geneInfo, 1, function(data){
					   if(data[7] == 4){
						   return(names(which(data == F)))
					   }else{
						   return(NA)
					   }
				 })) %>% mutate(Unique = replace(Unique, is.na(Unique), ""), Absent = replace(Absent, is.na(Absent), ""), Unique = paste0(Unique,Absent)) %>%
		select(Gene, ClustersPresent, Unique) %>% mutate(ClustersPresent = ifelse(ClustersPresent == 4, "Absent", "Present")) %>% 
		rename(PresentAbsent = ClustersPresent, AccessoryCluster = Unique)

geneInfo %>% group_by(AccessoryCluster) %>% count(PresentAbsent) %>% pivot_wider(id_cols = AccessoryCluster, names_from = PresentAbsent, values_from = `n`, values_fill = 0) %>%
	xtable(auto = T) %>% print(file = "UniquePresentAbsentTable.tex")

pdf("SharedAccessorybyCluster.pdf", width = 9, height = 9)
plot(euler(sharedGenesbyCluster, shape = "ellipse"), quantities = T,  legend = list(side = "bottom", nrow = 3, ncol = 2), fills = list(fill = clusterColours[names(sharedGenesbyCluster)]))
dev.off()

tmp <- coord %>% left_join(metaData %>% mutate(ST = as.character(ST)), by = c("Genome" = "Sample")) %>%
       	group_by(clusters) %>% count(Country, name = "Genomes") %>% 
	mutate(Prop = Genomes/sum(Genomes), Country = replace(Country, is.na(Country), "Unknown"), Region = groupedCountries[Country]) %>%
	mutate(Region = factor(Region, levels = names(countryColours)))

tmpLabels <- tmp %>% group_by(clusters) %>% summarize(Genome = sum(Genomes))

ggplot(tmp, aes(x = clusters,y = Prop, fill = Region)) +
       	geom_bar(stat = "identity") +
       	geom_text(data = tmpLabels,inherit.aes = F, aes(label = Genome, y = 1.025, x = clusters)) +
	annotate("rect", xmax = 0.55, xmin = 1.45, ymax = 0.5, ymin = 0.217, lty = 2, colour = "green", fill = NA, lwd = 1.5) +
	scale_fill_manual(values = unlist(countryColours)) +
	scale_y_continuous(breaks = pretty_breaks(n = 10)) +
	xlab("") + ylab("Proportion of Genomes") +
	theme(axis.text.x = element_text(angle = 11.75, hjust = 1, vjust = 1))
ggsave("CountryProportionsClustersGrouped.pdf", width = 6, height = 4)

######## Trying out a Corelation plot ###########
#accessDataSubDF <- as_tibble(accessDataSub)
#accessDataSubDF$Genome <- rownames(accessDataSub)
#
#occurenceDF <- accessDataSubDF %>% pivot_longer(-Genome, names_to = "Gene", values_to = "Count") %>%
#	left_join(metaData %>% filter(complete.cases(ST)), by  = c("Genome" = "Sample")) %>%
#	left_join(coord[,c("Genome", "clusters")]) %>%
#	mutate(Country = replace(Country, is.na(Country), "Unknown"), Region = groupedCountries[Country]) %>%
#	mutate(Region = factor(Region, levels = names(countryColours))) %>% select(-Country) %>% filter(Count != 0) %>%
#	group_by(clusters, Gene) %>% summarize(Count = sum(Count)) %>% pivot_wider(values_from = Count, names_from = Gene, values_fill = 0)
#
#clusterCount <- accessDataSubDF %>% filter(Genome != "Brancorsini") %>%
#	left_join(coord[,c("Genome", "clusters")]) %>% 
#	group_by(clusters) %>% count(name = "ClusterSize")
#
#clusterProp <- occurenceDF %>% pivot_longer(-clusters, names_to = "Gene", values_to = "Count") %>% 
#	left_join(coord %>% group_by(clusters) %>% count(name = "ClusterSize")) %>%
#	distinct() %>%
##	group_by(clusters, Gene) %>% mutate(ClusterProportion = Count/ClusterSize)
#	group_by(clusters, Gene) %>% summarize(ClusterProportion = Count/ClusterSize)
#
#annotationHeat <- clusterProp %>% mutate(ClusterProportion = ClusterProportion > 0.5, ClusterProportion = ifelse(ClusterProportion, "Present", "Absent")) %>% 
#       	pivot_wider(values_from = ClusterProportion, names_from = clusters) %>% as.data.frame()
#rownames(annotationHeat) <- annotationHeat$Gene
#annotationHeat <- annotationHeat[,-1]
#
#corMatrixColours <- list("Africa/America" = c("Present" = unname(clusterColours[3]), "Absent" = "white"), "Fertile Crescent" = c("Present" = unname(clusterColours[2]), "Absent" = "white"),
#			 "Indo-Pacific" = c("Present" = unname(clusterColours[4]), "Absent" = "white"), "Russia" = c("Present" = unname(clusterColours[5]), "Absent" = "white"),
#			 "Western Mediterranean" = c("Present" = unname(clusterColours[1]), "Absent" = "white"))
#
#occurenceMatrix <- occurenceDF[,-1] %>% as.matrix()
#rownames(occurenceMatrix) <- occurenceDF$clusters
#
#corMatrix <- cor(occurenceMatrix)
#
#pheatmap(corMatrix, clustering_method = "ward.D2",
#	 show_rownames = F, show_colnames = F, legend = T,
#	 cutree_rows = 5, cutree_cols = 5,
#	 annotation_row = annotationHeat, annotation_col = annotationHeat,
#	 annotation_legend = T, annotation_colors = corMatrixColours)
#
#treeCorMatrix <- hclust(as.dist(1 - corMatrix))
#plot(color_branches(treeCorMatrix, k = 5))

rm(clusterColoursTmp, NIPHUniqueClust3, notWmed, p1Clusts, p2clutss, tmp, tmp2, tmp3, uniqueAfAmGenes, tmpColours, ind, maxCount, coord, contrib, colour, coreColours, clustered, clusteredResults, annotationHeat, ann_colorsHeat, afAm, accessDataSub, accessData)

################################################
### Let's take a look at the Virulence Genes ###
################################################
# These are all coming from the Głowacka et al. 2018 paper Brucella - Virulence Factors, Pathogenesis and Treatment
identifiedCOGGenes <- read.delim("../BlastResults/Pangenome.tab") %>% as_tibble()

# Let's look for the genes, first in the names than the gene name itself
fromCOGName <- grep("LPS|lipopolisaccharide|T4SS|Type IV|VirB|VjbR|sodA|sodC|superoxide|sod|opg|urease|ure|ahpC|ahpD|ahp|nord|Nitric Oxide Reductase|Alkyl hydroperoxide reductase|cbb|cytochrome oxidase|Brucella virulence factor A|BvfA|xthA|Exonuclease III|Base excision repair|bvr|ompR",ignore.case = T, identifiedCOGGenes$Name)
fromCOGName <- identifiedCOGGenes$Query[fromCOGName]

# From Gene Name
fromGeneName <- grep("ompR|bvrr|bvrs|xth|bvf|nord|ahp|cco|ure|opg|sod|vjbr|virb|lps", ignore.case = T, roaryOutput$Gene)
fromGeneName <- roaryOutput$Gene[fromGeneName]

#
virulenceDf <- roaryOutput %>% filter(Gene %in% unique(c(fromCOGName,fromGeneName, "group_1825"))) # <--- group_1825 is bvfA per a blast run against only that gene

# Preparing the data
virulenceData <- as.data.frame(virulenceDf)
rownames(virulenceData) <- virulenceData[,1]
virulenceData <- virulenceData[,-1]
genomes <- colnames(virulenceData)
virulenceData <- apply(virulenceData, MARGIN = 1,FUN = as.numeric)
rownames(virulenceData) <- genomes

# We'll be removing the genes that are too common
virulenceDataSub <- virulenceData
maxCount <- virulenceDataSub %>% colSums() %>% max()
ind <- virulenceDataSub %>% colSums() >= (maxCount - 5)
virulenceDataSub <- virulenceDataSub[,!ind]

# Saving the genes that were shared amongst them ALL except for Geridu
tmp <- virulenceDf %>% as.data.frame()
rownames(tmp) <- tmp[,1]
tmp <- tmp[,-1]
genomes <- colnames(tmp)
tmp <- apply(tmp, MARGIN = 1,FUN = as.numeric)
rownames(tmp) <- genomes
maxCount <- tmp %>% colSums() %>% max()
colnames(tmp)[which(colSums(tmp) == 323)] %>% write.table("SharedVirGenes.list", col.names = F, row.names =F, quote = F)
######################################
# If I want to use a MDS plot
geneDist <- dist(virulenceDataSub, method = "binary")

fit <- cmdscale(geneDist, eig = T, k = 4, add =T)
coord <- fit$points %>% as_tibble()
coord$Genome <- rownames(fit$points)
contrib <- fit$eig/sum(fit$eig) * 100
coord <- coord %>% left_join(MLSTResults)
coord$ST[which(coord$Genome == "Brancorsini")] <- "Brancorsini"

# Clustering based on PCoA Coordinates
#fviz_nbclust(coord[,1:4],FUNcluster = clara, method = "silhouette")
#fviz_nbclust(coord[,1:4],FUNcluster = clara, method = "wss")
#
#clustered <- clara(coord[,-c(5,6)], 5)

#coord$clusters <- factor(clustered$clustering)
#write.table(coord, file = "FullPhyloClustering.tab", sep = "\t", row.names = F, col.names = T)

ann_colors$ST <- ann_colors$ST[mixedsort(unique(coord$ST))]
p1 <- coord %>%
       	ggplot(aes(x = V1, y = V2, label = Genome, colour = ST)) +
	geom_hline(yintercept=0, lty = 2, colour = "grey90") + 
	geom_vline(xintercept=0, lty = 2, colour = "grey90") +
	geom_point() +
	#stat_ellipse(show.legend = F, colour = "black") +
	xlab(bquote("PCoA 1 ("~.(round(contrib[1],2))~"%)")) +
	ylab(bquote("PCoA 2 ("~.(round(contrib[2],2))~"%)")) +
	scale_colour_manual(values = ann_colors$ST) +
	guides(colour = guide_legend(nrow = 2)) +
	#geom_text_repel(show.legend = F) +
	theme_classic() +
	theme(legend.position = "bottom")

p1

#############################################
# Saving the list of virulence genes found in the ancient genomes
virAncientOnly <- roaryOutput %>% filter(Gene %in% unique(c(fromCOGName, fromGeneName,"group_1825")), Brancorsini == 1) %>% select(Gene, Brancorsini)
virAncientOnly %>% write.table(file = "VirulenceGenes.tab", sep = "\t", row.names = F, quote = F)

#ancientOnly %>% filter(Gene %in% unique(c(fromCOGName, fromGeneName,"group_2139"))) %>% write.table(file = "VirulenceGenes.tab", sep = "\t", row.names =F, quote = F) # <-- 2139 comes from a blast run against bvfA
virAncientOnly <- virAncientOnly %>% left_join(identifiedCOGGenes, by = c("Gene" = "Query")) %>% rename(COGGene = Gene.y)
lipoPolyName <- virAncientOnly %>% filter(grepl("LPS|lipopolisaccharide", ignore.case =T, Name)| grepl("lps", ignore.case = T, COGGene) | grepl("lps", ignore.case = T, Gene)) %>%
	summarize(Genes = length(Gene)) %>% mutate(Category = "Lipopoly")
t4ss <- virAncientOnly %>% filter(grepl("T4SS|Type IV|VirB|VjbR", ignore.case =T, Name) | grepl("vjbr|virb", ignore.case = T, COGGene)| grepl("vjbr|virb", ignore.case = T, Gene)) %>%
	summarize(Genes = length(Gene)) %>% mutate(Category = "T4SS")
sod <- virAncientOnly %>% filter(grepl("sodA|sodC|superoxide|sod", ignore.case =T, Name) | grepl("sod", ignore.case = T, COGGene)| grepl("sod", ignore.case = T, Gene)) %>%
	summarize(Genes = length(Gene)) %>% mutate(Category = "SuperOxide")
cycicOPG <- virAncientOnly %>% filter(grepl("opg", ignore.case =T, Name) | grepl("opg", ignore.case = T, COGGene) | grepl("opg", ignore.case = T, Gene)) %>%
	summarize(Genes = length(Gene)) %>% mutate(Category = "OPG")
urease <- virAncientOnly %>% filter(grepl("urease|ure", ignore.case =T, Name) | grepl("ure", ignore.case = T, COGGene) | grepl("ure", ignore.case = T, Gene)) %>%
	summarize(Genes = length(Gene)) %>% mutate(Category = "Urease")
cytoOxi <- virAncientOnly %>% filter(grepl("cbb|cytochrome oxidase", ignore.case =T, Name)| grepl("cco", ignore.case = T, COGGene) | grepl("cco", ignore.case = T, Gene)) %>%
	summarize(Genes = length(Gene)) %>% mutate(Category = "CytoOxi")
ahp <- virAncientOnly %>% filter(grepl("ahpC|ahpD|ahp|Alkyl hydroperoxide reductase", ignore.case =T, Name) | grepl("ahp", ignore.case = T, COGGene) | grepl("ahp", ignore.case = T, Gene)) %>%
	summarize(Genes = length(Gene)) %>% mutate(Category = "Alkyl")
nor <- virAncientOnly %>% filter(grepl("nord|Nitric Oxide Reductase", ignore.case =T, Name) | grepl("nord", ignore.case = T, COGGene) | grepl("nord", ignore.case = T, Gene)) %>%
	summarize(Genes = length(Gene)) %>% mutate(Category = "Nitric")
bruceVir <- virAncientOnly %>% filter(grepl("Brucella virulence factor A|BvfA", ignore.case =T, Name) | grepl("bvf|group_1825", ignore.case = T, COGGene) | grepl("bvf|group_1825", ignore.case = T, Gene)) %>%
	summarize(Genes = length(Gene)) %>% mutate(Category = "Brucella Virulence Factor")
exo <- virAncientOnly %>% filter(grepl("xthA|Exonuclease III|Base excision repair", ignore.case =T, Name) | grepl("xth", ignore.case = T, COGGene) | grepl("xth", ignore.case = T, Gene)) %>%
	summarize(Genes = length(Gene)) %>% mutate(Category = "Exonuclease")
bvr <- virAncientOnly %>% filter(grepl("bvr|ompR", ignore.case =T, Name) | grepl("ompR|bvrr|bvrs|bvf", ignore.case = T, COGGene) | grepl("ompR|bvrr|bvrs|bvf", ignore.case = T, Gene)) %>%
	summarize(Genes = length(Gene)) %>% mutate(Category = "BvrR/BvrS")

virCatSummarized <- bind_rows(list(lipoPolyName, t4ss, sod, cycicOPG, urease, cytoOxi, ahp, nor,bruceVir, exo, bvr))
xtable(as.data.frame(virCatSummarized)) %>% print(file = "VirulenceGenesSummarized.tex", include.rownames = F)
write.csv(file = "VirulenceGenesSummarize.csv", virCatSummarized, row.names = F, quote = F)
write.table(unique(c(fromCOGName, fromGeneName)), file = "VirulenceGenes.list", quote = F, row.names = F, col.names = F)

###################################
### Now to search the AMR Genes ###
###################################

cardBlastnucl <- BlastParsing("AMRData/NuclearBlast.tab", 7) %>% rename(Gene = Query) %>%
       	separate(Match, into = c("DB", "NCBI Accession", "Strand", "Location", "ARO", "AROGene"), sep = "\\|") %>% mutate(ARO = gsub("ARO:","", ARO))
cardAROinfo <- fromJSON(file = "card.json")[1:4967] # Removing the metadata data at the end

foundAMRinfo <- pblapply(cardAROinfo, cl = 7,function(x){ # This pulls out the actual ARO info that we need. Can delete cardAROinfo afterwards
	if(x$ARO_accession %in% cardBlastnucl$ARO){
		return(x)
	}
			 })
rm(cardAROinfo)
foundAMRinfo <- foundAMRinfo[!sapply(foundAMRinfo, is.null)] # Removing all the nulls....
names(foundAMRinfo) <- unname(sapply(foundAMRinfo, function(x) return(x$ARO_accession))) # Renaming the elements so that I can pull them out via ARO accession IDs

foundAncientAMR <- ancientOnly %>% filter(Gene %in% cardBlastnucl$Gene)
amrCat <- lapply(foundAMRinfo, function(x){ # Pulling out the ARO_cats
	       lapply(x$ARO_category, function(y){
		      return(y$category_aro_name)
		})
			 })
amrCatNames <- names(amrCat)

for (x in amrCatNames){amrCat[[x]] <- amrCat[[x]] %>% unlist() %>% unname()} # Making it a single list...
rm(amrCatNames)

cardBlastnucl <- cardBlastnucl %>% select(Gene, ARO, AROGene, PIdent, Evalue) %>% 
	mutate(AROCategory = amrCat[ARO] %>% unname())

ancientAMR <- ancientOnly %>% inner_join(cardBlastnucl) %>% select(Sample, Gene, MeanCoverage, Name, AROGene)
xtable(ancientAMR) %>% print(file = "AMRTable.tex", include.rownames = F)

################
# Finally, looking at the core genes which were absent. What's their COG Functions?
################
cogDefinitions <- read.delim("../BlastResults/fun-20.tab", header =F , col.names = c("ID", "Colour", "Function")) # Requires a download of the COG 2020 DB
cogDefList <- as.list(cogDefinitions %>% pull(Function))
names(cogDefList) <- cogDefinitions$ID

branFunctionalCategories <- identifiedCOGGenes %>% filter(Query %in% coreGenesMissingBrancorsini$Gene) %>% pull(FunctionalCategory) %>% strsplit(split = "") %>% unlist() %>% table() %>% as_tibble()
colnames(branFunctionalCategories) <- c("Function", "Brancorsini")
branFunctionalCategories <- branFunctionalCategories %>% mutate(Function = unlist(cogDefList[Function])) %>% arrange(-Brancorsini)

branFunctionalCategories %>% xtable() %>% print(file = "MissingCore.tex", include.rownames =F)

