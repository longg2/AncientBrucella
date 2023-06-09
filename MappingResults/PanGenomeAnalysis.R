library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(ape)
library(scales)
library(ggExtra)
library(ggnewscale)
library(cluster)
library(parallel)
library(pbapply)
library(purrr)
library(pheatmap)
library(factoextra)
library(gtools)
library(xtable)
library(eulerr)
library(rjson)
##### Functions #####
VCFParsing <- function(vcf){
		# Preparing for the final part
		FinalName <- gsub(".*\\/|\\.vcf", "",vcf)
	# Getting the file Ready
	vcfFile <- read.delim(file = vcf,comment.char = "#", header = F,
			      			      col.names = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "unknown")) %>%
				filter(FILTER == "PASS") %>% as_tibble()
					#filter(REF != "C" & ALT != "T") %>% filter(REF != "G" & ALT != "A") %>% as <- tibble() # Filtered out the G -> A and C -> T Transitions
			
				if(nrow(vcfFile) == 0){
							cat("No SNPs found in",FinalName, "\n")
		return(NA)
			}
					
					tmp <- unlist(strsplit(vcfFile$FORMAT[1], ":"))
						tmp2 <- unlist(strsplit(gsub("=.*?(?=;)|=.*$","",vcfFile$INFO[1], perl = T), ";"))
			vcfFile <- vcfFile %>% select(-FORMAT) %>% separate(unknown, into = tmp, sep = ":") %>%
						mutate(INFO = gsub("(?<=;).*?=|^.*?=", "", INFO, perl = T)) %>% separate(INFO, into= tmp2, sep = ";")
						vcfFile[,c(6,8:11,14)]<- sapply(vcfFile[,c(6,8:11,14)], as.numeric)
							return(vcfFile)
}

CorrectingGC <- function(depthData){ # From the iREP Paper
	GCModel <- lm(Mean ~ GCContent, data = depthData)
	resids <- abs(summary(GCModel)$residuals) # Trying to find the top 1% largest residuals. Don't care about location
	residFilt <- resids < quantile(resids,0.99)
	GCModelFilt <- lm(Mean ~ GCContent, data = depthData[residFilt,])
	adjr2 <- summary(GCModel)$adj.r.squared
	beforeGC <- depthData[residFilt,] %>% ggplot(aes(x = GCContent, y = Mean)) +
		geom_point() + geom_density_2d(colour = "grey") + theme_classic() + geom_smooth(method = "lm") +
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
		geom_point() + geom_density_2d(colour = "grey") + theme_classic() + geom_smooth(method = "lm") +
	       	ggtitle(bquote("After GC Correction"~R[adj]^2 == .(round(adjr2, 3))), subtitle = "Top 1% Residuals Filtered") +
	       	ylab("Mean Coverage") + xlab("GC Content")
	GCPlots <- ggarrange(beforeGC, afterGC, ncol = 2, align = "hv")
	return(list("Plot" = GCPlots, UpdatedDepths = depthData))

}

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
##### Some Basic Setup #####
theme_set(theme_classic())

colour <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6',
'#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3',
'#808000', '#ffd8b1', '#000075', '#808080', '#f0f0f0', '#000000')

clusterColours <- c("Western Mediterranean" = "#4D9DE0", "Fertile Crescent" = "#E15554", "Africa/America" = "#E1BC29", "Indo-Pacific" = "#3BB273", "Russia" = "#7768AE")
tmp <- read.delim("CountriesGrouped.list", header = F, col.names = c("Country", "Colour"))
#tmp <- read.delim("CountriesInterestAlphabetized.list", header = F, col.names = c("Country", "Colour"))
countryColours <- as.list(tmp$Colour)
names(countryColours) <- tmp$Country
rm(tmp)

# Reading the MLST Results
MLSTResults <- read.delim("MLSTResults.txt") %>% as_tibble()
MLSTResults <- MLSTResults %>% select(Sample, ST) %>% mutate(ST = gsub("\\*| ","", ST))
colnames(MLSTResults)[1] <- "Genome"

# Reading the metadata table
metaData <- read.delim("MetadataAll.tab", header =T) %>% as_tibble
ann_colors <- list(ST = c("5" = colour[3], "7" = colour[15], "8" = colour[10], "9" = colour[9], "10" = colour[6], "11" = colour[17], "12" = colour[8], "39" = colour[13], "40" = colour[2], "41" = colour[16], "42" = colour[7], "43" = colour[4], "71" = colour[11], "88" = colour[12], "102" = colour[19],"Ancient" = colour[1], "NF" = colour[20], "NIPH" = colour[22]),
	Norway = c("Normal" = "#00205B" , "Odd" = "#BA0C2F", "Other" = "#FFFFFF"))

############ Let's get the Depths ###########
op <- pboptions(type = "timer")
ncores = 8
files <- list.files("Depths", full.names = T)
depthDf <- do.call(bind_rows, pblapply(files, cl = ncores,function(f){
		tmp <- as_tibble(read.delim(f, header = F, col.names = c("Gene", "Pos", "Coverage")))
		tmp$Genome <- gsub(".*/", "", gsub("\\..*|_genomic","",f)) 	
		tmp <- tmp %>% group_by(Genome, Gene) %>%
			summarize(MeanCoverage = mean(Coverage), sdCoverage = sd(Coverage), PercentCoverage = sum(Coverage > 0)/length(Coverage), .groups = "drop") %>%
			mutate(CV = sdCoverage/MeanCoverage)
		return(tmp)
}))


# Getting the GC Content of the Genes
gcContent <- read.table("PanGenomeGC.tab", header = F, col.names = c("Gene", "Length","GC"), sep = "\t") %>% as_tibble() %>% mutate(GC = GC/100)
depthDf <- depthDf %>% left_join(gcContent)

# Now to correct for GC Bias
GCBiased <- depthDf %>% filter(MeanCoverage > 0) %>% select(Genome,Gene,MeanCoverage, GC)
colnames(GCBiased) <- c("Genome","Gene","Mean", "GCContent")
GCCorrected <- lapply(split(GCBiased, GCBiased$Genome), function(x)CorrectingGC(x)$UpdatedDepths) %>% bind_rows()
depthDf <- depthDf %>% full_join(GCCorrected %>% select(Genome,Gene, GCCorrected)) %>% mutate(GCCorrected = replace(GCCorrected, is.na(GCCorrected), 0))

######### Figuring out the thresholds ##############
CVSearch <- lapply(seq(0,5,0.5), function(x) depthDf %>% filter(GCCorrected > 0, CV <= x) %>% count(Genome, name = as.character(x))) %>%
       	reduce(full_join, by = "Genome") %>% pivot_longer(-Genome, names_to = "CV", values_to = "Genes") %>%
       	mutate(Genes = replace(Genes, is.na(Genes), 0),DiffGenes = diff(c(0, Genes)), CV = as.numeric(CV)) %>%
	filter(CV > 0)

CVResults <- CVSearch %>% ggplot(aes(x = CV, y = DiffGenes)) +
	theme_classic() +
	geom_line() + geom_point() + ylab("New Genes") +
	geom_vline(xintercept = 1.5, lty = 2, colour = colour[1]) +
	xlab("CV Threshold") 
#ggsave(CVResults, "CVSearch.pdf", width = 6, height = 4)

HistRect <- depthDf %>% filter(CV <= 1.5, GCCorrected > 0) %>% summarize(Mean = mean(GCCorrected), SD = sd(GCCorrected))

histPlotBoth <- depthDf %>% filter(CV <= 1.5, GCCorrected > 0) %>%
	ggplot(aes(x = GCCorrected)) +
#	geom_rect(data = HistRect, inherit.aes = F,
#		  aes(xmin = ifelse(Mean - 2 * SD >= 0,Mean - 2 * SD,0) , xmax = Mean + 2 * SD, ymin = -Inf, ymax = Inf),
#		  colour = colour[4], lty = 1, fill = NA) +
#	geom_rect(data = HistRect, inherit.aes = F,
#		  aes(xmin = ifelse(Mean - 3 * SD >= 0,Mean - 3 * SD,0) , xmax = Mean + 3 * SD, ymin = -Inf, ymax = Inf),
#		  colour = colour[5], lty = 1, fill = NA) +
	geom_histogram(position = "identity", colour = "black", fill = colour[1], binwidth = 1) +
	geom_vline(xintercept = HistRect$Mean + c(-2 *HistRect$SD, -3 * HistRect$SD), colour = colour[4:5]) +
	#geom_vline(xintercept = HistRect$Mean + c(2 *HistRect$SD, -2 * HistRect$SD), colour = colour[4]) +
	#geom_vline(xintercept = HistRect$Mean + c(3 *HistRect$SD, -3 * HistRect$SD), colour = colour[5]) +
	geom_vline(data = HistRect, aes(xintercept = Mean), colour = "black", lty = 2) +
	xlab("Mean Read Depth") + ylab("Genes") + theme(legend.position = "bottom") + 
	scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) + scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) 
ggarrange(CVResults, histPlotBoth, labels = "auto", align = "hv")
ggsave("CVSearch.pdf", width = 6, height = 4)

threshBoth <- lapply(seq(0,5,0.5), function(x){
		lapply(seq(0,20,0.5), function(y){
				tmp <- depthDf %>% filter(CV <= x, GCCorrected >= y) %>% count(Genome, name = "Count") %>% pull(Count)
				df <- tibble(CV = x, Mean = y, Count = tmp)
				return(df)
			}) %>% bind_rows()
		}) %>% bind_rows()

threshBoth %>% ggplot(aes(x = CV, y = Mean, fill = Count)) +
	geom_tile(colour = "black") +
	geom_text(aes(x = CV, y = Mean, label = Count), colour = ifelse(threshBoth %>% pull(Count) >= 2996, "black", "black")) +
	geom_rect(aes(xmin = 1.25, xmax = 5.25, ymin = -0.25, ymax = 7.25), fill = NA, colour = colour[1], lty = 1, lwd = 1) +
	geom_rect(aes(xmin = 1.25, xmax = 1.75, ymin = 11.25, ymax = 11.75), fill = NA, colour = colour[4], lty = 1, lwd = 1) +
	geom_rect(aes(xmin = 1.25, xmax = 1.75, ymin = 6.75, ymax = 7.25), fill = NA, colour = colour[5], lty = 1, lwd = 1) +
	scale_fill_gradient2(low = colour[4], mid = "white", high = colour[1],
			     midpoint = 2996,
			     #midpoint = min(CVSearch$Count) + (max(CVSearch$Count) - min(CVSearch$Count))/2,
			     breaks = breaks_pretty(n = 5),
			     limits = c(min(threshBoth$Count), max(threshBoth$Count))) +
	scale_x_continuous(breaks = breaks_pretty(n = 5)) +
	scale_y_continuous(breaks = breaks_pretty(n = 10)) +
	theme(legend.position = "bottom")

ggsave("ThresholdTest.pdf", width = 12, height = 9)

############# Core Gene Presence #####
roaryOutput <- as_tibble(read.delim("PresenceAbsence.tab"))
colnames(roaryOutput)[2:ncol(roaryOutput)] <- gsub("^X|\\.scaffold|\\.genome|\\.result|_genomic", "", colnames(roaryOutput)[2:ncol(roaryOutput)])
colnames(roaryOutput)[2:ncol(roaryOutput)] <- gsub("\\.", "-", colnames(roaryOutput)[2:ncol(roaryOutput)])
# Finding the Core Genome
coreGenes <- rowSums(roaryOutput[,-1]) >= floor(0.99 * ncol(roaryOutput))
coreGenes <- roaryOutput$Gene[coreGenes]

########### Some basic Coverage Data  ###############
depthDf %>% mutate(Status = ifelse(Gene %in% coreGenes, "Core", "Accessory")) %>% filter(Genome == "Brancorsini") %>% filter(CV <= 1.5) %>%
       	mutate(Status = factor(Status, levels = c("Core","Accessory"))) %>%
	ggplot(aes(x = GCCorrected, fill = Status)) + 
	geom_histogram(colour = "black", binwidth = 1) +
	geom_vline(xintercept = HistRect$Mean + c(-2 *HistRect$SD, -3 * HistRect$SD), colour = colour[4:5]) +
	geom_vline(data = HistRect, aes(xintercept = Mean), colour = "black", lty = 2) +
	#geom_vline(xintercept = CorsiniRect$Mean, colour = "purple", lty = 2) +
	scale_fill_manual(values = c(Accessory = colour[4],Core = colour[1]), "Genome") +
	#scale_fill_manual(values = c(Accessory = "#007dba",Core = "#1b998b")) +
	xlab("Mean Read Depth") + ylab("Genes") + theme(legend.position = "bottom") + 
	scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) + scale_y_continuous(breaks = scales::pretty_breaks(n = 10))
	#theme(axis.title.x = element_blank(), axis.text.x = element_blank())

# Core genes found in ancient sample
SD2 <- HistRect$Mean - 2 * HistRect$SD
SD3 <- HistRect$Mean - 3 * HistRect$SD
genesAncient <- depthDf %>% filter(CV <= 1.5, GCCorrected >= SD2) %>% pull(Gene)
length(coreGenes) - sum(coreGenes %in% genesAncient)

genesAncient <- depthDf %>% filter(CV <= 1.5, GCCorrected >= SD3) %>% pull(Gene)
length(coreGenes) - sum(coreGenes %in% genesAncient)

depthDf %>% filter(Gene %in% coreGenes[!(coreGenes %in% genesAncient)])

threshBoth <- lapply(seq(0,5,0.5), function(x){
		lapply(seq(0,20,0.5), function(y){
				tmp <- depthDf %>% filter(CV <= x, GCCorrected >= y, Gene %in% coreGenes) %>% count(Genome, name = "Count") %>% pull(Count)
				df <- tibble(CV = x, Mean = y, Count = tmp)
				return(df)
			}) %>% bind_rows()
		}) %>% bind_rows()

threshBoth %>% ggplot(aes(x = CV, y = Mean, fill = Count)) +
	geom_tile(colour = "black") +
	geom_text(aes(x = CV, y = Mean, label = Count)) +
	geom_rect(aes(xmin = 1.75, xmax = 5.25, ymin = -0.25, ymax = 3.75), fill = NA, colour = colour[1], lty = 1, lwd = 1) +
	geom_rect(aes(xmin = 1.25, xmax = 1.75, ymin = 11.25, ymax = 11.75), fill = NA, colour = colour[4], lty = 1, lwd = 1) +
	geom_rect(aes(xmin = 1.25, xmax = 1.75, ymin = 6.75, ymax = 7.25), fill = NA, colour = colour[5], lty = 1, lwd = 1) +
	scale_fill_gradient2(low = colour[4], mid = "white", high = colour[1],
			     midpoint = 2826,
			     #midpoint = min(CVSearch$Count) + (max(CVSearch$Count) - min(CVSearch$Count))/2,
			     breaks = breaks_pretty(n = 3),
			     limits = c(min(threshBoth$Count), max(threshBoth$Count))) +
	scale_x_continuous(breaks = breaks_pretty(n = 5)) +
	scale_y_continuous(breaks = breaks_pretty(n = 10)) +
	theme(legend.position = "bottom")

ggsave("ThresholdTestCore.pdf", width = 12, height = 9)

# Some quick t-tests
tmp <- depthDf %>% filter(Genome == "Brancorsini") %>% mutate(Status = ifelse(Gene %in% coreGenes, "Core", "Accessory")) %>% 
	filter(CV <= 1.5, GCCorrected >= SD2)
t.test(tmp %>% filter(Status == "Core") %>% pull(GCCorrected), tmp %>% filter(Status != "Core") %>% pull(GCCorrected))

tmp <- depthDf %>% filter(Genome == "Brancorsini") %>% mutate(Status = ifelse(Gene %in% coreGenes, "Core", "Accessory")) %>% 
	filter(CV <= 1.5, GCCorrected >= SD3)
t.test(tmp %>% filter(Status == "Core") %>% pull(GCCorrected), tmp %>% filter(Status != "Core") %>% pull(GCCorrected))

##### Gene Presence table #####
SD2Genes <- depthDf %>% filter(Genome == "Brancorsini", CV <= 1.5, GCCorrected >= SD2) %>% select(-c(CV,MeanCoverage,sdCoverage, PercentCoverage, Genome)) %>%
	mutate(.keep = "unused", AncientSD2 = GCCorrected)

SD3Genes <- depthDf %>% filter(Genome == "Brancorsini", CV <= 1.5, GCCorrected >= SD3) %>% select(-c(CV,MeanCoverage,sdCoverage, PercentCoverage, Genome)) %>%
	mutate(.keep = "unused", AncientSD3 = GCCorrected)

AllPA <- roaryOutput %>% left_join(SD2Genes %>% select(-c(Length, GC))) %>% left_join(SD3Genes %>% select(-c(Length, GC))) %>%
       	mutate(AncientSD2 = ifelse(is.na(AncientSD2)|AncientSD2 == 0, 0,1)) %>%
       	mutate(AncientSD3 = ifelse(is.na(AncientSD3)|AncientSD3 == 0, 0,1))

# Merging the Corsini and Geridu Data
ancientOnly <- SD2Genes %>% 
	pivot_longer(cols = AncientSD2, names_to = "Sample", values_to = "MeanCoverage", values_drop_na = T) %>%
	mutate(Status = ifelse(Gene %in% coreGenes, "Core", "Accessory")) %>%
       	inner_join(depthDf %>% select(Gene, GC) %>% distinct())

ancientOnly %>% 
	 summarize(Length = sum(Length), GCMean = mean(GC), GCSE = qnorm(0.975) * sd(GC)/sqrt(length(GC)), GCHi = GCMean + GCSE, GCLo = GCMean - GCSE, MeanDepth = mean(MeanCoverage), MeanSE = qnorm(0.975)*sd(MeanCoverage)/sqrt(length(MeanCoverage)), MeanHi = MeanDepth + MeanSE, MeanLo = MeanDepth - MeanSE) %>% as.data.frame()
 ancientOnly %>% count(Status)

tmp <- AllPA %>% select(-Gene) %>% summarize_all(sum) %>% t() %>% as.data.frame()
tmp$Genome <- rownames(tmp)
geneCounts <- tmp %>% as_tibble() 
geneCounts <- geneCounts %>% left_join(MLSTResults)
colnames(geneCounts)[1] <- "Genes"
geneCounts %>% filter(!grepl("Ancient", Genome)) %>% ggplot(aes(y = Genes, x ="", color = ST)) +
       	geom_boxplot(colour = "black") + 
	#geom_jitter(height = 0) +
	geom_point(data = geneCounts %>% filter(grepl("Ancient", Genome)), aes(shape = Genome, x = ""), colour = colour[1]) +
	scale_colour_manual(values = colour[4:5]) + 
	theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), legend.position = "right") +
	guides(colour = guide_legend(title = "Sample"))
ggsave("GeneCountsBmelFixed.png", width = 6, height = 9)

##### Core Gene Scatter Plot #####
depthDf %>% filter(Gene %in% coreGenes) %>%
	ggplot(aes(y = CV, x = GCCorrected)) + geom_point() +
	geom_vline(xintercept = SD3, lty =2 , colour = "red") +
	geom_hline(yintercept = 1.5, lty = 2, colour = "red") 

# Going to do this in a more sane way
coreDf <- AllPA %>% filter(Gene %in% coreGenes)

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
coord$ST[(nrow(coord)-2):nrow(coord)] <- "Brancorsini"

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
ggsave(pCore, file = "PCoACoreSept2022.pdf", width = 6, height = 4)

##### Accessory Genome Analysis #####
accessDf <- AllPA %>% filter(!(Gene %in% coreGenes))

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

# If I want to use a MDS plot
geneDist <- dist(accessDataSub, method = "binary")

fit <- cmdscale(geneDist, eig = T, k = 4, add =T)
coord <- fit$points %>% as_tibble()
coord$Genome <- rownames(fit$points)
contrib <- fit$eig/sum(fit$eig) * 100
coord <- coord %>% left_join(MLSTResults)
coord$ST[grepl("Ancient", coord$Genome)] <- c("Ancient")
coord$Shape <- ifelse(grepl("Ancient.*", coord$Genome), ifelse(grepl("AncientSD2", coord$Genome), "2SD", "3SD"), "Modern")
#coord$ST[coord$Genome %in% c("Brancorsini", "Geridu")] <- c("Brancorsini", "Geridu")

# Clustering based on PCoA Coordinates
fviz_nbclust(coord[,1:4],FUNcluster = clara, method = "silhouette")
fviz_nbclust(coord[,1:4],FUNcluster = clara, method = "wss")
#fviz_nbclust(coord[,1:4],FUNcluster = clara, method = "gap_stat")
#
clustered <- clara(coord[,-c(5,6)], 5)

coord$clusters <- clustered$clustering
coord  <- coord %>% mutate(clusters = ifelse(clusters == 1, "Western Mediterranean",
				 ifelse(clusters == 2, "Fertile Crescent",
				       	ifelse(clusters == 3, "Africa/America", ifelse(clusters == 4, "Indo-Pacific", "Russia")))))
#write.table(coord, file = "FullPhyloClustering.tab", sep = "\t", row.names = F, col.names = T)
coreColours <- ann_colors$ST[names(ann_colors$ST) %in% accessPlotData$ST]
coreColours <- coreColours[mixedsort(names(coreColours))]

accessPlotData <- coord %>% left_join(metaData %>% select(-ST), by = c("Genome" = "Sample")) %>%
       	#mutate(Norway3 = ifelse(grepl("^NIPH-*|NI_2007", Genome) & grepl("8",ST), "OddNorway", NA)) %>%
       	mutate(ST = replace(ST, grepl("^NIPH-*|NI_2007", Genome) & grepl("8",ST), "NIPH")) %>%
       	mutate(ST = factor(ST, levels = names(coreColours)))

# Getting the labels for the ellipses
labelCoord <- accessPlotData %>% group_by(clusters) %>% summarize(V1 = mean(V1), V2 = mean(V2), V3 = mean(V3), V4 = mean(V4))

p1 <- accessPlotData %>%
       	ggplot(aes(x = V1, y = V2, label = Genome, colour = ST, group = clusters, shape = Shape)) +#, group = clusters)) +
	geom_hline(yintercept=0, lty = 2, colour = "grey90") + 
	geom_vline(xintercept=0, lty = 2, colour = "grey90") +
	geom_point() +
	scale_shape_manual(values = c("Modern" = 16, "2SD" = 15, "3SD" = 17)) +
	scale_colour_manual(values = coreColours, "Sequence Type") +
	guides(colour = guide_legend(nrow = 3, title.position = "top", title.hjust = 0.5)) +
	new_scale_colour() +
	stat_ellipse(mapping = aes(colour = clusters)) +
	scale_colour_manual("Accessory Gene Clustering", values = clusterColours) +
	#stat_ellipse(data = . %>% filter(!is.na(Norway3)),show.legend = F, colour = "black", lty = 2,mapping = aes(group = Norway3)) +
	xlab(bquote("PCoA 1 ("~.(round(contrib[1],2))~"%)")) +
	ylab(bquote("PCoA 2 ("~.(round(contrib[2],2))~"%)")) +
	geom_text_repel(data = labelCoord, aes(x = V1, y = V2, label = clusters, colour = clusters), inherit.aes = F, show.legend = F) +
	guides(colour = guide_legend(nrow = 3, title.position = "top", title.hjust = 0.5),
	       shape = guide_legend(nrow = 2, title.position = "top", title.hjust = 0.5))  +
	theme(legend.position = "bottom")
	#theme(legend.position = "bottom", axis.text.x = element_blank(), axis.title.x = element_blank())

p1
#STLabelledPCoA <- ggarrange(p1,p2, nrow = 1, align = "hv", common.legend = T, legend = "bottom", labels = "auto")
ggsave(file = "ReviewMappingPCoA.pdf", width = 9, height = 6)

# Now to look at a heatmap
annotationHeat <- accessPlotData %>% select(Genome, ST, clusters) %>% distinct() %>% as.data.frame()
rownames(annotationHeat) <- annotationHeat$Genome
annotationHeat <- annotationHeat[,-1]
colnames(annotationHeat)[2] <- "Accessory Clusters"

ann_colorsHeat <- list("Accessory Clusters" = clusterColours, ST = ann_colors$ST)

pheatmap(accessDataSub, clustering_distance_rows = "binary",clustering_distance_cols = "binary", clustering_method = "ward.D2",
	 show_rownames = F, show_colnames = F, legend = F,
	 cutree_rows = 5,
	 annotation_row = annotationHeat, annotation_colors = ann_colorsHeat, annotation_legend = T,
	 annotation_names_row = F, filename = "AccessoryHeatmapPresentation.pdf", width = 6, height = 6)

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
pdf("SharedAccessorybyCluster.pdf", width = 9, height = 9)
plot(euler(sharedGenesbyCluster, shape = "ellipse"), quantities = T,  legend = list(side = "bottom", nrow = 3, ncol = 2), fills = list(fill = clusterColours[names(sharedGenesbyCluster)]))
dev.off()
#ggvenn(sharedGenesbyCluster, fill_color = c("#E1BC29", "#E15554", "#3BB273", "#7768AE", "#4D9DE0"))
#ggsave("SharedAccessorybyCluster.pdf", width = 9, height = 9)

groupedCountries <- c("Afghanistan" = "Asia", "Argentina" = "Americas", "China" = "China", "Ethiopia" = "Africa", "Georgia" = "Asia", "India" = "India", "Iraq" = "Fertile Crescent","Israel" = "Fertile Crescent", "Italy" = "Italy", "Kuwait" = "Fertile Crescent", "Malta" = "Europe", "Nigeria" = "Africa", "Norway" = "Norway", "Protugal" = "Europe","Somalia" = "Africa", "Turkey" = "Fertile Crescent", "United Kingdom" = "Europe", "Zimbabwe" = "Africa", "Unknown" = "Unknown", "Albania" = "Europe","Bulgaria" = "Europe", "Canada" = "Americas", "Cyprus" = "Europe", "Iran" = "Fertile Crescent", "Jordan" = "Fertile Crescent", "Kosovo" = "Europe","Malaysia" = "Asia", "Pakistan" = "Asia", "Russia" = "Russia", "Saudi Arabia" = "Asia", "Sudan" = "Africa", "Syria" = "Fertile Crescent", "Thailand" = "Asia","Turkmenistan" = "Asia", "Egypt" = "Egypt", "France" = "Europe", "Morocco" = "Morocco", "United States" = "Americas", "Portugal" = "Europe")

tmp <- coord %>% left_join(metaData %>% mutate(ST = as.character(ST)), by = c("Genome" = "Sample")) %>%
       	group_by(clusters) %>% count(Country, name = "Genomes") %>% 
	mutate(Prop = Genomes/sum(Genomes), Country = replace(Country, is.na(Country), "Unknown"), Region = groupedCountries[Country]) %>%
	mutate(Region = factor(Region, levels = names(countryColours)))

tmpLabels <- tmp %>% group_by(clusters) %>% summarize(Genome = sum(Genomes))

ggplot(tmp, aes(x = clusters,y = Prop, fill = Region)) +
       	geom_bar(stat = "identity") +
       	geom_text(data = tmpLabels,inherit.aes = F, aes(label = Genome, y = 1.025, x = clusters)) +
	scale_fill_manual(values = unlist(countryColours)) +
	scale_y_continuous(breaks = pretty_breaks(n = 10)) +
	xlab("") + ylab("Proportion of Genomes") +
	theme(axis.text.x = element_text(angle = 11.75, hjust = 1, vjust = 1))
ggsave("CountryProportionsClustersGrouped.pdf", width = 6, height = 4)

#ggplot(tmp, aes(x = "",y = Prop, fill = Country)) + geom_bar(stat = "identity", width = 1) + coord_polar("y", start = 0) + facet_wrap("clusters")

tmp <- coord %>% left_join(metaData %>% mutate(ST = as.character(ST)), by = c("Genome" = "Sample")) %>%
       	group_by(clusters) %>% count(Host, name = "Genomes") %>% 
	mutate(Prop = Genomes/sum(Genomes), Host = replace(Host, is.na(Host), "Unknown")) %>%
	mutate(clusters = ifelse(clusters == 1, "Western Mediterranean",
				 ifelse(clusters == 2, "Fertile Crescent",
				       	ifelse(clusters == 3, "Africa/America", ifelse(clusters == 4, "Indo-Pacific", "Russia")))))
ggplot(tmp, aes(x = clusters,y = Prop, fill = Host)) +
       	geom_bar(stat = "identity") +
       	geom_text(data = tmpLabels,inherit.aes = F, aes(label = Genome, y = 1.025, x = clusters)) +
	#scale_fill_manual(values = unlist(countryColours)) +
	scale_y_continuous(breaks = pretty_breaks(n = 10)) +
	xlab("") + ylab("Proportion of Genomes") +
	theme(axis.text.x = element_text(angle = 22.5, hjust = 1, vjust = 1))

coord %>% select(Genome, clusters) %>% write.table(file = "GenomeAccessClustered.tab", sep = "\t", row.names = F)

rm(clusterColoursTmp, NIPHUniqueClust3, notWmed, p1Clusts, p2clutss, tmp, tmp2, tmp3, uniqueAfAmGenes, tmpColours, ind, maxCount, coord, contrib, colour, coreColours, clustered, clusteredResults, annotationHeat, ann_colorsHeat, afAm, accessDataSub, accessData)

##### Let's take a look at the Virulence Genes #####
# These are all coming from the Głowacka et al. 2018 paper Brucella - Virulence Factors, Pathogenesis and Treatment
identifiedCOGGenes <- read.delim("../BlastResults/COGClassified90/Pangenome.tab") %>% as_tibble()

# Let's look for the genes, first in the names than the gene name itself
fromCOGName <- grep("LPS|lipopolisaccharide|T4SS|Type IV|VirB|VjbR|sodA|sodC|superoxide|sod|opg|urease|ure|ahpC|ahpD|ahp|nord|Nitric Oxide Reductase|Alkyl hydroperoxide reductase|cbb|cytochrome oxidase|Brucella virulence factor A|BvfA|xthA|Exonuclease III|Base excision repair|bvr|ompR",ignore.case = T, identifiedCOGGenes$Name)
fromCOGName <- identifiedCOGGenes$Query[fromCOGName]

# From Gene Name
fromGeneName <- grep("ompR|bvrr|bvrs|xth|bvf|nord|ahp|cco|ure|opg|sod|vjbr|virb|lps", ignore.case = T, AllPA$Gene)
fromGeneName <- AllPA$Gene[fromGeneName]

#
virulenceDf <- AllPA %>% filter(Gene %in% unique(c(fromCOGName,fromGeneName, "group_2139"))) # <--- group_2139 is bvfA per a blast run against only that gene

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
tmp <- virulenceDf[,-which(colnames(virulenceDf) == "Geridu")] %>% as.data.frame()
rownames(tmp) <- tmp[,1]
tmp <- tmp[,-1]
genomes <- colnames(tmp)
tmp <- apply(tmp, MARGIN = 1,FUN = as.numeric)
rownames(tmp) <- genomes
maxCount <- tmp %>% colSums() %>% max()
colnames(tmp)[which(colSums(tmp) == 324)] %>% write.table("SharedVirGenes.list", col.names = F, row.names =F, quote = F)
######################################
# If I want to use a MDS plot
geneDist <- dist(virulenceDataSub, method = "binary")

fit <- cmdscale(geneDist, eig = T, k = 4, add =T)
coord <- fit$points %>% as_tibble()
coord$Genome <- rownames(fit$points)
contrib <- fit$eig/sum(fit$eig) * 100
coord <- coord %>% left_join(MLSTResults)
coord$ST[(nrow(coord)-1):nrow(coord)] <- c("Brancorsini", "Geridu")

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

##### Now to do to the SNP heterozygosity analysis #####
ancientOnly %>% group_by(Sample) %>% summarize(sum(Length)) # What's our estimated genome lengths?
ancientOnly <- ancientOnly %>% group_by(Sample) %>% mutate(CopyNumber = MeanCoverage/mean(MeanCoverage))
geneLengths <- gcContent %>% select(-GC)

vcf <- VCFParsing("HetData/Brancorsini.vcf.gz") %>% select(CHROM,FILTER, QUAL, GT)

# Getting the total gene set
vcfHet <- vcf %>% mutate(Hetero = ifelse(grepl("0/1|1/0", GT), T, F))

corsiniGenes <- geneLengths %>% filter(Gene %in% Corsini$Gene) # Gets the relevant Genes
vcfCorsini <- vcfHet %>% filter(!is.na(FILTER),QUAL >= 100) %>% # does all of the heavy lifting
       	group_by(CHROM) %>% summarize(Hetero = sum(Hetero)) %>% right_join(corsiniGenes, by = c("CHROM" = "Gene")) %>%
	mutate(Hetero = replace(Hetero, is.na(Hetero),0)) %>% group_by(CHROM) %>% summarize(HeteroFrac = Hetero/Length) %>%
	mutate(HeteroFrac = replace(HeteroFrac, is.infinite(HeteroFrac), 0), Status = ifelse(CHROM %in% coreGenes, "Core", "Accessory")) %>%
	mutate(Status = factor(Status, levels = c("Core", "Accessory")))%>% mutate(Sample = "Brancorsini")

coreHet <-vcfCorsini %>% filter(Status == "Core") %>% pull(HeteroFrac)
accessHet <-vcfCorsini %>% filter(Status == "Accessory") %>% pull(HeteroFrac)
car:::Anova(type = 3, lm(HeteroFrac ~ Status, vcfCorsini))

## Now for Geridu
#vcf <- VCFParsing("HetData/Nodule1.vcf.gz") %>% select(CHROM,FILTER, QUAL, GT)
#
#vcfHet <- vcf %>% mutate(Hetero = ifelse(grepl("0/1|1/0", GT), T, F))
#
#geriduGenes <- geneLengths %>% filter(Gene %in% Geridu$Gene)
#vcfGeridu <- vcfHet %>% filter(!is.na(FILTER),QUAL >= 100) %>%
#       	group_by(CHROM) %>% summarize(Hetero = sum(Hetero)) %>% right_join(corsiniGenes, by = c("CHROM" = "Gene")) %>%
#	mutate(Hetero = replace(Hetero, is.na(Hetero),0)) %>% group_by(CHROM) %>% summarize(HeteroFrac = Hetero/Length) %>%
#	mutate(HeteroFrac = replace(HeteroFrac, is.infinite(HeteroFrac), 0), Status = ifelse(CHROM %in% coreGenes, "Core", "Accessory")) %>%
#	mutate(Status = factor(Status, levels = c("Core", "Accessory"))) %>% mutate(Sample = "Geridu")
#
#coreHet <-vcfGeridu %>% filter(Status == "Core") %>% pull(HeteroFrac)
#accessHet <-vcfGeridu %>% filter(Status == "Accessory") %>% pull(HeteroFrac)
#car:::Anova(type = 3, lm(HeteroFrac ~ Status, vcfGeridu))
#
# Summary Statistics of both
vcfPlot <- vcfCorsini# %>% bind_rows(vcfGeridu)
vcfPlot %>% group_by(Sample) %>%
       	summarize(Mean = mean(HeteroFrac), SD = sd(HeteroFrac), confInt = qnorm(0.975)*SD/sqrt(length(HeteroFrac)), 
			 high = Mean + confInt, low = Mean - confInt) %>% as.data.frame()
vcfPlot %>% group_by(Sample,Status) %>%
       	summarize(Mean = mean(HeteroFrac), SD = sd(HeteroFrac), confInt = qnorm(0.975)*SD/sqrt(length(HeteroFrac)), 
			 high = Mean + confInt, low = Mean - confInt) %>% as.data.frame()

coreHetNo <- vcfPlot %>% filter(Sample == "Brancorsini") %>% filter(HeteroFrac == 0, Status == "Core") %>%
	summarize(Genes =length(CHROM))
accessHetNo <- vcfPlot  %>% filter(Sample == "Brancorsini") %>% filter(HeteroFrac == 0, Status == "Accessory") %>%
	summarize(Genes =length(CHROM))

violinPoints <- vcfPlot %>% filter(Sample == "Brancorsini") %>% group_by(Status) %>%
       	summarize(Mean = mean(HeteroFrac), SD = sd(HeteroFrac), confInt = qnorm(0.975)*SD/sqrt(length(HeteroFrac)), 
			 high = Mean + confInt, low = Mean - confInt) %>% as.data.frame()
violinCounts <- vcfPlot %>% filter(Sample == "Brancorsini") %>% group_by(Status) %>%
	summarize(Genes = length(Status), GenesWithHet = sum(HeteroFrac > 0))

hetHist <- vcfPlot %>% filter(Sample == "Brancorsini") %>% ggplot(aes(x = Status, y= HeteroFrac, fill = Status)) +
       #	geom_vline(xintercept = vcfPlot %>% filter(HeteroFrac < 0) %>% summarize(mean(HeteroFrac)) %>% pull(), colour = "red", lty = 2) +
	geom_violin(draw_quantiles = 0.5) +
	geom_errorbar(data = violinPoints, inherit.aes = F, mapping = aes(ymin = low, ymax = high, x = Status), width = 0.15) +
	geom_point(data = violinPoints, inherit.aes = F, mapping = aes(y = Mean, x = Status)) +
	geom_text(violinCounts, inherit.aes = F, mapping = aes(y = 10^-4, x = Status, label = Genes)) +
	geom_text(violinCounts, inherit.aes = F, mapping = aes(y = 3*10^-5, x = Status, label = GenesWithHet), colour = "red") +
	scale_y_log10() + annotation_logticks(sides = "l") +
	ylab("P(Heterozygous)") + xlab("") +
	scale_fill_manual(values = c(Accessory = colour[4],Core = colour[1])) +
	theme(legend.position = "bottom") 
#hetHist <- vcfPlot %>% filter(Sample == "Brancorsini") %>% ggplot(aes(x = HeteroFrac, fill = Status)) +
#       	geom_histogram(position = "identity", alpha = 0.75, colour = "black") +
#       #	geom_vline(xintercept = vcfPlot %>% filter(HeteroFrac < 0) %>% summarize(mean(HeteroFrac)) %>% pull(), colour = "red", lty = 2) +
#       	theme_classic() +
#	scale_x_log10() + annotation_logticks(sides = "b") + scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
#	ylab("Genes") + xlab("P(Heterozygous)") +
#	scale_fill_manual(values = c(Core = "#f8333c", Accessory = "#007dba")) +
#	theme(legend.position = "bottom") +
#	annotate(geom = "text", label = paste0("No Heterozygous SNPs: ",coreHetNo$Genes), x = 3e-3, y = 14, colour = "#f8333c") +
#	annotate(geom = "text", label = paste0("No Heterozygous SNPs: ",accessHetNo$Genes), x = 3e-3, y = 12, colour = "#007dba")

# Now let's get the copy numbers involved
ancientOnly <- ancientOnly %>% left_join(vcfPlot %>% select(-Status), by = c("Sample", "Gene" = "CHROM"))
       	
model <- lm(data = ancientOnly %>% filter(HeteroFrac > 0, Sample == "Brancorsini"), log10(HeteroFrac) ~ CopyNumber*Status)
r2 <- round(summary(model)$adj.r.squared,2)
tmp <- summary(model)$fstatistic
pval <- round(pf(tmp[1],tmp[2],tmp[3], lower.tail = F),3)

ancientOnly %>% filter(Sample == "Brancorsini") %>% ggplot(aes(x = GC, y = HeteroFrac, colour = Status)) +
	geom_point() +
	theme_bw() +
	geom_smooth(method = "lm") +
	scale_y_log10() + annotation_logticks(sides = "l") + scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
	ylab("P(Heterozygous)") +
	scale_colour_manual(values = c(Accessory = colour[4],Core = colour[1])) +
	theme(legend.position = "bottom") 

copyHet <- ancientOnly %>% filter(Sample == "Brancorsini") %>% ggplot(aes(x = CopyNumber, y = HeteroFrac, colour = Status)) +
	geom_point() +
	theme_bw() +
	geom_smooth(method = "lm") +
	scale_y_log10() + annotation_logticks(sides = "l") + scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
	xlab("Copy Number") + ylab("P(Heterozygous)") +
	scale_colour_manual(values = c(Accessory = colour[4],Core = colour[1])) +
	annotate(geom = "text", x = 1.6, y = 3e-3, label = bquote(R[adj]^2 == .(r2))) +
	annotate(geom = "text", x = 1.6, y = 2e-3, label = bquote(P %~~% .(pval))) +
	theme(legend.position = "bottom") 

write.table(ancientOnly %>% filter(HeteroFrac > 0), file = "HeterozygousGenes.tab", sep = "\t", quote = F, row.names = F )
ggsave(copyHet, file = "CopyNumberHeterozygosity.png", width = 6, height = 4)

hetPlots <- ggarrange(hetHist, copyHet, ncol = 1, common.legend = T, legend = "bottom", align = "hv", labels = "AUTO")

tmp <- ggplot() + theme_void()

left <- ggarrange(GCPlot, nrow = 1, labels = c("a"), align = "hv")
right <- ggarrange(histPlotBran, hetHist, common.legend = T, legend = "bottom", nrow = 1, labels = c("b","c"), align = "hv")

ggarrange(left, right, ncol = 1)
ggsave("Fig2Sept.pdf", width = 9, height = 6)
#ggarrange(histPlotBoth,histPlotBran, tmp, hetHist, common.legend = T, legend = "bottom", labels = "auto")


ggarrange(tmp, hetPlots, common.legend = T)
ggsave(file = "HeterozygosityJune2022.pdf", height = 6, width = 9)

hetHist %>% ggsave(file = "~/Documents/University/LabMeetings/2022/Mar25/Figures/HetHist.pdf", width = 6, height= 4)
copyHet %>% ggsave(file = "~/Documents/University/LabMeetings/2022/Mar25/Figures/HetCN.pdf", width = 6, height= 4)

# Looking for potential HGTs and Transposons
ancientOnly <- ancientOnly %>% left_join(read.delim("../BlastResults/COGClassified90/Pangenome.tab"), by = c("Gene" = "Query"))
ancientOnly %>% filter(Sample == "Brancorsini") %>% filter(HeteroFrac > 0, Status == "Accessory") %>% select(Gene, CopyNumber,HeteroFrac, Name)

# Saving the high copy number genes
highCopyCor <- ancientOnly %>% filter(Sample == "Brancorsini", MeanCoverage > (mean(MeanCoverage) + 2*sd(MeanCoverage)))
highCopyGer <- ancientOnly %>% filter(Sample != "Brancorsini", MeanCoverage > (mean(MeanCoverage) + 2*sd(MeanCoverage)))
highCopyCor %>% bind_rows(highCopyGer) %>% write.table(file = "HighCopyGenes.tab", sep = "\t", row.names =F, quote = F)
#highCopyCor %>% bind_rows(highCopyGer) %>% pull(Gene) %>% unique() %>% write.table(file = "HighCopyNumberGenes.list", sep = "\t", row.names =F, quote = F)

# 
tmp <- ancientOnly %>% filter(Sample == "Brancorsini") %>% mutate(High = ifelse(MeanCoverage > (mean(MeanCoverage) + 2*sd(MeanCoverage)), T, F))
t.test(tmp %>% filter(High, HeteroFrac > 0) %>% pull(HeteroFrac), tmp %>% filter(!High) %>% pull(HeteroFrac))

#############################################
# Saving the list of virulence genes found in the ancient genomes
ancientOnly %>% filter(Gene %in% unique(c(fromCOGName, fromGeneName,"group_2139"))) %>% write.table(file = "VirulenceGenes.tab", sep = "\t", row.names =F, quote = F) # <-- 2139 comes from a blast run against bvfA

virAncientOnly <- ancientOnly %>% filter(Gene %in% unique(c(fromCOGName, fromGeneName, "group_2139"))) 

lipoPolyName <- virAncientOnly %>% filter(grepl("LPS|lipopolisaccharide", ignore.case =T, Name) | grepl("lps", ignore.case = T, Gene)) %>%
	summarize(Mean = mean(MeanCoverage), Genes = length(Gene), Error = qnorm(0.975) * sd(MeanCoverage)/sqrt(Genes), ErrorLow = Mean - Error, ErrorHigh = Mean + Error) %>% mutate(Category = "Lipopoly")
t4ss <- virAncientOnly %>% filter(grepl("T4SS|Type IV|VirB|VjbR", ignore.case =T, Name) | grepl("vjbr|virb", ignore.case = T, Gene)) %>%
	summarize(Mean = mean(MeanCoverage), Genes = length(Gene), Error = qnorm(0.975) * sd(MeanCoverage)/sqrt(Genes), ErrorLow = Mean - Error, ErrorHigh = Mean + Error) %>% mutate(Category = "T4SS")
sod <- virAncientOnly %>% filter(grepl("sodA|sodC|superoxide|sod", ignore.case =T, Name) | grepl("sod", ignore.case = T, Gene)) %>%
	summarize(Mean = mean(MeanCoverage), Genes = length(Gene), Error = qnorm(0.975) * sd(MeanCoverage)/sqrt(Genes), ErrorLow = Mean - Error, ErrorHigh = Mean + Error) %>% mutate(Category = "SuperOxide")
cycicOPG <- virAncientOnly %>% filter(grepl("opg", ignore.case =T, Name) | grepl("opg", ignore.case = T, Gene)) %>%
	summarize(Mean = mean(MeanCoverage), Genes = length(Gene), Error = qnorm(0.975) * sd(MeanCoverage)/sqrt(Genes), ErrorLow = Mean - Error, ErrorHigh = Mean + Error) %>% mutate(Category = "OPG")
urease <- virAncientOnly %>% filter(grepl("urease|ure", ignore.case =T, Name) | grepl("ure", ignore.case = T, Gene)) %>%
	summarize(Mean = mean(MeanCoverage), Genes = length(Gene), Error = qnorm(0.975) * sd(MeanCoverage)/sqrt(Genes), ErrorLow = Mean - Error, ErrorHigh = Mean + Error) %>% mutate(Category = "Urease")
cytoOxi <- virAncientOnly %>% filter(grepl("cbb|cytochrome oxidase", ignore.case =T, Name) | grepl("cco", ignore.case = T, Gene)) %>%
	summarize(Mean = mean(MeanCoverage), Genes = length(Gene), Error = qnorm(0.975) * sd(MeanCoverage)/sqrt(Genes), ErrorLow = Mean - Error, ErrorHigh = Mean + Error) %>% mutate(Category = "CytoOxi")
ahp <- virAncientOnly %>% filter(grepl("ahpC|ahpD|ahp|Alkyl hydroperoxide reductase", ignore.case =T, Name) | grepl("ahp", ignore.case = T, Gene)) %>%
	summarize(Mean = mean(MeanCoverage), Genes = length(Gene), Error = qnorm(0.975) * sd(MeanCoverage)/sqrt(Genes), ErrorLow = Mean - Error, ErrorHigh = Mean + Error) %>% mutate(Category = "Alkyl")
nor <- virAncientOnly %>% filter(grepl("nord|Nitric Oxide Reductase", ignore.case =T, Name) | grepl("nord", ignore.case = T, Gene)) %>%
	summarize(Mean = mean(MeanCoverage), Genes = length(Gene), Error = qnorm(0.975) * sd(MeanCoverage)/sqrt(Genes), ErrorLow = Mean - Error, ErrorHigh = Mean + Error) %>% mutate(Category = "Nitric")
bruceVir <- virAncientOnly %>% filter(grepl("Brucella virulence factor A|BvfA", ignore.case =T, Name) | grepl("bvf|group_2139", ignore.case = T, Gene)) %>%
	summarize(Mean = mean(MeanCoverage), Genes = length(Gene), Error = qnorm(0.975) * sd(MeanCoverage)/sqrt(Genes), ErrorLow = Mean - Error, ErrorHigh = Mean + Error) %>% mutate(Category = "Brucella Virulence Factor")
exo <- virAncientOnly %>% filter(grepl("xthA|Exonuclease III|Base excision repair", ignore.case =T, Name) | grepl("xth", ignore.case = T, Gene)) %>%
	summarize(Mean = mean(MeanCoverage), Genes = length(Gene), Error = qnorm(0.975) * sd(MeanCoverage)/sqrt(Genes), ErrorLow = Mean - Error, ErrorHigh = Mean + Error) %>% mutate(Category = "Exonuclease")
bvr <- virAncientOnly %>% filter(grepl("bvr|ompR", ignore.case =T, Name) | grepl("ompR|bvrr|bvrs|bvf", ignore.case = T, Gene)) %>%
	summarize(Mean = mean(MeanCoverage), Genes = length(Gene), Error = qnorm(0.975) * sd(MeanCoverage)/sqrt(Genes), ErrorLow = Mean - Error, ErrorHigh = Mean + Error) %>% mutate(Category = "BvrR/BvrS")

virCatSummarized <- bind_rows(list(lipoPolyName, t4ss, sod, cycicOPG, urease, cytoOxi, ahp, nor,bruceVir, exo, bvr))
xtable(as.data.frame(virCatSummarized)) %>% print(file = "VirulenceGenesSummarized.tex", include.rownames = F)
#virCatSummarized %>% xtable(auto = T) %>% print(file = "~/Documents/University/ComprehensiveExamination/Paper/VirTable.tex")
write.csv(file = "VirulenceGenesSummarize.csv", virCatSummarized, row.names = F, quote = F)
write.table(unique(c(fromCOGName, fromGeneName)), file = "../VirulenceGenes.list", quote = F, row.names = F, col.names = F)

##### Now to search the AMR Genes #####

cardBlastnucl <- BlastParsing("AMRData/NuclearBlast.tab", 7) %>% rename(Gene = Query) %>%
       	separate(Match, into = c("DB", "NCBI Accession", "Strand", "Location", "ARO", "AROGene"), sep = "\\|") %>% mutate(ARO = gsub("ARO:","", ARO))
cardAROinfo <- fromJSON(file = "../card.json")[1:4967] # Removing the metadata data at the end

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

##### Finally, looking at the core genes which were absent. What's their COG Functions?####
cogDefinitions <- read.delim("../BlastResults/fun-20.tab", header =F , col.names = c("ID", "Colour", "Function"))
cogDefList <- as.list(cogDefinitions %>% pull(Function))
names(cogDefList) <- cogDefinitions$ID

branFunctionalCategories <- identifiedCOGGenes %>% filter(Query %in% branMissingCore) %>% pull(FunctionalCategory) %>% strsplit(split = "") %>% unlist() %>% table() %>% as_tibble()
colnames(branFunctionalCategories) <- c("Function", "Brancorsini")
branFunctionalCategories <- branFunctionalCategories %>% mutate(Function = unlist(cogDefList[Function]))

gerFunctionalCategories <- identifiedCOGGenes %>% filter(Query %in% gerMissingCore) %>% pull(FunctionalCategory) %>% strsplit(split = "") %>% unlist() %>% table() %>% as_tibble()
colnames(gerFunctionalCategories) <- c("Function", "Geridu")
gerFunctionalCategories %>% mutate(Function = unlist(cogDefList[Function])) %>% full_join(branFunctionalCategories) %>% 
	arrange(-Brancorsini, -Geridu) %>% select(Function, Brancorsini, Geridu) %>% xtable() %>% print(file = "MissingCore.tex", include.rownames =F)

