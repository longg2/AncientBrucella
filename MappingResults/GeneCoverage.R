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

############ Now to actually plot some of the missing genes #########
SD2 <- HistRect$Mean - 2 * HistRect$SD
SD3 <- HistRect$Mean - 3 * HistRect$SD
SD4 <- HistRect$Mean - 4 * HistRect$SD

genesAncient <- depthDf %>% filter(CV <= 1.5, GCCorrected >= SD3) %>% pull(Gene)
missingCore <- depthDf %>% filter(Gene %in% coreGenes[!(coreGenes %in% genesAncient)])

# Missing Genes
tmp <- as_tibble(read.delim(files[1], header = F, col.names = c("Gene", "Pos", "Coverage"))) %>% filter(Gene %in% coreGenes)

tmp %>% filter(Gene %in% c("truB","group_2323", "group_1907")) %>% 
	ggplot(aes(x = Pos, y = Coverage)) +
	geom_line(colour = colour[1]) +
	scale_y_continuous(trans = "log2")  +
	theme_classic() +
	facet_wrap("Gene", scales = "free_x")
ggsave("MissingCoreGeneDepths.pdf", width = 6, height = 4)
