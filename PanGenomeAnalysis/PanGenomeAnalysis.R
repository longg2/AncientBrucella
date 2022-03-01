library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)
#library(ggvenn)
library(scales)
library(ggExtra)
#library(reshape2)
library(cluster)
library(parallel)
library(pbapply)
library(purrr)
#library(pheatmap)
#library(dendextend)
#library(FactoMineR)
#library(factoextra)
library(gtools)

## Functions
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
MLSTResults <- MLSTResults %>% select(Sample, ST) %>% mutate(ST = gsub("\\*| ","", ST))
colnames(MLSTResults)[1] <- "Genome"

#ann_colors <- list(ST = c("9" = colour[9], "11" = colour[17], "43" = colour[4], "88" = colour[12], "Geridu" = colour[5],"Corsini" = colour[1], "NF" = colour[20], "Reference" = colour[22]))
ann_colors <- list(ST = c("5" = colour[10], "7" = colour[15], "8" = colour[3], "9" = colour[9], "10" = colour[6], "11" = colour[17], "12" = colour[8], "39" = colour[13], "40" = colour[2], "41" = colour[16], "42" = colour[7], "43" = colour[4], "71" = colour[11], "88" = colour[12], "102" = colour[19], "Geridu" = colour[5],"Brancorsini" = colour[1], "NF" = colour[20], "Reference" = colour[22]))

############ Let's get the Depths ###########
op <- pboptions(type = "timer")
ncores = 8
depthDf <- do.call(bind_rows, pblapply(c("Depths/JessSamples.tab.gz", "Depths/KayBMel.tab.gz"), cl = ncores,function(f){
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
depthDf <- depthDf %>% mutate(Genome = gsub("Trimmed","",Genome))

############# Core Gene Presence #####
roaryOutput <- as_tibble(read.delim("PresenceAbsence.tab"))
colnames(roaryOutput)[2:ncol(roaryOutput)] <- gsub("X|\\.scaffold|\\.genome|\\.result|_genomic", "", colnames(roaryOutput)[2:ncol(roaryOutput)])
colnames(roaryOutput)[2:ncol(roaryOutput)] <- gsub("\\.", "-", colnames(roaryOutput)[2:ncol(roaryOutput)])
# Finding the Core Genome
coreGenes <- rowSums(roaryOutput[,-1]) >= floor(0.95 * ncol(roaryOutput))
#coreGenes <- rowSums(ecoliOnly[,-1]) >= floor(0.99 * ncol(ecoliOnly))
coreGenes <- roaryOutput$Gene[coreGenes]

########### Some basic Coverage Data  ###############
Corsini <- depthDf %>% filter(Genome == "JessSamples", CV <=1) %>%
	pivot_wider(names_from = c(Genome), values_from=GCCorrected) %>% select(-c(MeanCoverage,GC,sdCoverage, PercentCoverage, CV)) #%>%
	#mutate(Genome = "10x - CV Filter")
CorsiniRect <- Corsini %>% summarize(Mean = mean(JessSamples), SD2 = 2*sd(JessSamples))
Geridu <- depthDf %>% filter(Genome != "JessSamples", CV <=1) %>%
	pivot_wider(names_from = c(Genome), values_from=GCCorrected) %>% select(-c(MeanCoverage,GC,sdCoverage, PercentCoverage, CV))

GeriduRect <- Geridu %>% summarize(Mean = mean(KayBMel), SD2 = 2*sd(KayBMel))

corsiniPlot <- Corsini %>% mutate(Status = ifelse(Gene %in% coreGenes, "Core", "Accessory")) %>% ggplot(aes(x = JessSamples, fill = Status)) + 
	geom_rect(data = CorsiniRect, inherit.aes = F, aes(xmin = Mean - SD2, xmax = Mean + SD2, ymin = -Inf, ymax = Inf), colour = "purple", alpha = 0.5) +
	geom_histogram(position = "identity", alpha = 0.75, colour = "black", binwidth = 1) +
	geom_vline(xintercept = 10, colour = "black", lty = 2) +
	geom_vline(xintercept = CorsiniRect$Mean, colour = "purple", lty = 2) +
	scale_fill_manual(values = c(Accessory = "#007dba",Core = "#f8333c")) + theme_bw() +
	xlab("Mean Read Depth") + ylab("Genes") + theme(legend.position = "bottom", axis.title.x = element_blank(), axis.text.x = element_blank()) + 
	scale_x_continuous(breaks = scales::pretty_breaks(n = 10), limits = c(0,50)) + scale_y_continuous(breaks = scales::pretty_breaks(n = 10), limits = c(0, 1250))

geriduPlot <- Geridu %>% mutate(Status = ifelse(Gene %in% coreGenes, "Core", "Accessory")) %>% ggplot(aes(x = KayBMel, fill = Status)) + 
	geom_rect(data = GeriduRect, inherit.aes = F, aes(xmin = Mean - SD2, xmax = Mean + SD2, ymin = -Inf, ymax = Inf), colour = "purple", alpha = 0.5) +
	geom_histogram(position = "identity", alpha = 0.75, colour = "black", binwidth = 1) +
	geom_vline(xintercept = 1, colour = "black", lty = 2) +
	geom_vline(xintercept = GeriduRect$Mean, colour = "purple", lty = 2) +
	scale_fill_manual(values = c(Core = "#f8333c",Accessory = "#007dba")) + theme_bw() +
	xlab("Mean Read Depth") + ylab("Genes") + theme(legend.position = "bottom") + 
	scale_x_continuous(breaks = scales::pretty_breaks(n = 10), limits = c(0,50)) + scale_y_continuous(breaks = scales::pretty_breaks(n = 10), limits = c(0, 1250))

tmp <- ggarrange(corsiniPlot, geriduPlot, nrow = 2, labels= "AUTO", common.legend = T, align = "hv", legend = "bottom")

ggsave(tmp, file = "GenePresenceHistogram.pdf", width = 9, height = 6)

# Genes which were > +2SD from distribution
#####################################################################
# Gene Presence table
Corsini <- Corsini %>% filter(JessSamples >= 10)
Geridu <- Geridu %>% filter(KayBMel >= 1)
AllPA <- roaryOutput %>% left_join(full_join(Corsini, Geridu) %>% select(-Length)) %>%
       	mutate(JessSamples = ifelse(is.na(JessSamples)|JessSamples == 0, 0,1),
	KayBMel = ifelse(is.na(KayBMel)|KayBMel == 0, 0,1)) 
colnames(AllPA)[(ncol(AllPA)-1):ncol(AllPA)] <- c("Brancorsini", "Geridu")
colnames(AllPA) <- gsub("--","-XXX-", colnames(AllPA))

# Merging the Corsini and Geridu Data
ancientOnly <- Corsini %>% full_join(Geridu) %>% pivot_longer(cols = JessSamples:KayBMel, names_to = "Sample", values_to = "MeanCoverage", values_drop_na = T) %>%
	mutate(Sample = ifelse(Sample == "KayBMel", "Geridu", "Brancorsini"), Status = ifelse(Gene %in% coreGenes, "Core", "Accessory"))

# Now to make a boxplot to compare our gene to everyone else
tmp <- AllPA %>% select(-Gene) %>% summarize_all(sum) %>% t() %>% as.data.frame()
tmp$Genome <- rownames(tmp)
geneCounts <- tmp %>% as_tibble() 
geneCounts <- geneCounts %>% left_join(MLSTResults)
colnames(geneCounts)[1] <- "Genes"
geneCounts %>% filter(!grepl("Brancorsini|Geridu", Genome)) %>% ggplot(aes(y = Genes, x ="", color = ST)) +
       	geom_boxplot(outlier.shape = NA, colour = "black") + theme_bw() +
	#geom_jitter(height = 0) +
	geom_point(data = geneCounts %>% filter(grepl("Brancorsini|Geridu", Genome)), aes(colour = Genome, x = "")) +
	scale_colour_manual(values = ann_colors$ST) + 
	theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), legend.position = "right")
ggsave("GeneCountsBmel.pdf", width = 6, height = 9)
#####################################################################
# We want to compare the genes in the two samples
SamplePA <- Corsini %>% full_join(Geridu) %>% mutate(CorsiniPresent = ifelse(JessSamples >= 10, T, F), GeriduPresent = ifelse(KayBMel >= 1, T,F)) %>%
	mutate(CorsiniPresent = replace(CorsiniPresent, is.na(CorsiniPresent), F), GeriduPresent = replace(GeriduPresent, is.na(GeriduPresent), F)) %>%
	mutate(Status = ifelse(Gene %in% coreGenes, "Core", "Accessory"))

tmp <- coreGenes[!(coreGenes %in% SamplePA$Gene)]
tmp <- c(tmp,SamplePA %>% filter(Status == "Core", !CorsiniPresent, !GeriduPresent) %>% pull(Gene))
write.table(tmp,"CoreGenesMissingAll.list", col.names =F, row.names = F, quote =F)

SamplePA %>% filter(Status == "Accessory", CorsiniPresent, GeriduPresent == F) %>% pull(Gene) %>% write.table("AccessoryGenesMissingGeridu.list", col.names =F, row.names = F, quote =F)
SamplePA %>% filter(Status == "Accessory", CorsiniPresent == F, GeriduPresent) %>% pull(Gene) %>% write.table("AccessoryGenesMissingCorsini.list", col.names =F, row.names = F, quote =F)

table(SamplePA$CorsiniPresent, SamplePA$GeriduPresent)

#####################################################################
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
coord$ST[(ncol(coord)-1):ncol(coord)] <- c("Brancorsini", "Geridu")

# Quickly filtering the list so that only the STs present are used
pCore <- coord %>%
       	ggplot(aes(x = V1, y = V2, label = Genome, colour = ST)) +
	geom_hline(yintercept=0, lty = 2, colour = "grey90") + 
	geom_vline(xintercept=0, lty = 2, colour = "grey90") +
	geom_point() +
	xlab(bquote("PCoA 1 ("~.(round(contrib[1],2))~"%)")) +
	ylab(bquote("PCoA 2 ("~.(round(contrib[2],2))~"%)")) +
	scale_colour_manual(values = ann_colors$ST) +
	guides(colour = guide_legend(nrow = 2)) +
	#geom_text_repel(show.legend = F) +
	theme_classic() +
	theme(legend.position = "bottom") 
pCore
ggsave(pCore, file = "PCoACore.pdf", width = 6, height = 4)

########## Accessory Genome Analysis ########
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

######################################
# If I want to use a MDS plot
geneDist <- dist(accessDataSub, method = "binary")

fit <- cmdscale(geneDist, eig = T, k = 4, add =T)
coord <- fit$points %>% as_tibble()
coord$Genome <- rownames(fit$points)
contrib <- fit$eig/sum(fit$eig) * 100
coord <- coord %>% left_join(MLSTResults)
coord$ST[(ncol(coord)-1):ncol(coord)] <- c("Brancorsini", "Geridu")

# Clustering based on PCoA Coordinates
sil <- sapply(2:6, function(i){clara(coord[,-c(5,6)], i)$silinfo$avg.width})
plot(2:6, sil, type = "b")

clustered <- clara(coord[,-c(5,6)], 4)

coord$clusters <- factor(clustered$clustering)

ann_colors$ST <- ann_colors$ST[mixedsort(unique(coord$ST))]
p1 <- coord %>%
       	ggplot(aes(x = V1, y = V2, label = Genome, colour = ST, group = clusters)) +
	geom_hline(yintercept=0, lty = 2, colour = "grey90") + 
	geom_vline(xintercept=0, lty = 2, colour = "grey90") +
	geom_point() +
	#stat_ellipse(lty = 2,show.legend = F) +
	xlab(bquote("PCoA 1 ("~.(round(contrib[1],2))~"%)")) +
	ylab(bquote("PCoA 2 ("~.(round(contrib[2],2))~"%)")) +
	scale_colour_manual(values = ann_colors$ST) +
	guides(colour = guide_legend(nrow = 2)) +
	#geom_text_repel(show.legend = F) +
	theme_classic() +
	theme(legend.position = "bottom")

p1
ggsave(p1, file = "AccessoryPCoA.pdf", width = 9, height = 6)

# Pulling out the Genomes in cluster 1 and doing the same thing
italy <- coord %>% filter(clusters == 1) %>% pull(Genome)
geneDist <- dist(accessDataSub[rownames(accessDataSub) %in% italy,], method = "binary")

fit <- cmdscale(geneDist, eig = T, k = 4, add =T)
coord <- fit$points %>% as_tibble()
coord$Genome <- rownames(fit$points)
contrib <- fit$eig/sum(fit$eig) * 100
coord <- coord %>% left_join(MLSTResults)
coord$ST[(ncol(coord)-1):ncol(coord)] <- c("Brancorsini", "Geridu")

# Getting the clusters ready (if we want them)
sil <- sapply(2:6, function(i){clara(coord[,-c(5,6)], i)$silinfo$avg.width})
plot(2:6, sil, type = "b")

clustered <- clara(coord[,-c(5,6)], 3)
coord$clusters <- factor(clustered$clustering)

ann_colors$ST <- ann_colors$ST[mixedsort(unique(coord$ST))]
p1Clust <- coord %>%
       	ggplot(aes(x = V1, y = V2, label = Genome, colour = ST, group = clusters)) +
	geom_hline(yintercept=0, lty = 2, colour = "grey90") + 
	geom_vline(xintercept=0, lty = 2, colour = "grey90") +
	geom_point() +
	#stat_ellipse(show.legend = F) +
	xlab(bquote("PCoA 1 ("~.(round(contrib[1],2))~"%)")) +
	ylab(bquote("PCoA 2 ("~.(round(contrib[2],2))~"%)")) +
	scale_colour_manual(values = ann_colors$ST) +
	#guide(shape = )
	#geom_text_repel(show.legend = F) +
	theme_classic() +
	theme(legend.position = "bottom")
p1Clust
ggsave(p1Clust, file = "AccessoryItalyPCoA.pdf", width = 9, height = 6)

ggarrange(p1,p1Clust, legend = "bottom", align = "hv", nrow = 1, common.legend = T, labels = "AUTO")
ggsave(file = "PCoA_Accessory_ClusterWhole.pdf", width = 9, height = 6)

################################################
# Now to do to the SNP heterozygosity analysis #
################################################
ancientOnly %>% group_by(Sample) %>% summarize(sum(Length)) # What's our estimated genome lengths?
ancientOnly <- ancientOnly %>% group_by(Sample) %>% mutate(CopyNumber = MeanCoverage/mean(MeanCoverage))
geneLengths <- gcContent %>% select(-GC)

vcf <- VCFParsing("HetData/CorsiniHet.vcf") %>% select(CHROM,FILTER, QUAL, GT)

# Getting the total gene set
vcfHet <- vcf %>% mutate(Hetero = ifelse(grepl("0/1|1/0", GT), T, F))

corsiniGenes <- geneLengths %>% filter(Gene %in% Corsini$Gene) # Gets the relevant Genes
vcfCorsini <- vcfHet %>% filter(!is.na(FILTER),QUAL >= 30) %>% # does all of the heavy lifting
       	group_by(CHROM) %>% summarize(Hetero = sum(Hetero)) %>% right_join(corsiniGenes, by = c("CHROM" = "Gene")) %>%
	mutate(Hetero = replace(Hetero, is.na(Hetero),0)) %>% group_by(CHROM) %>% summarize(HeteroFrac = Hetero/Length) %>%
	mutate(HeteroFrac = replace(HeteroFrac, is.infinite(HeteroFrac), 0), Status = ifelse(CHROM %in% coreGenes, "Core", "Accessory")) %>%
	mutate(Status = factor(Status, levels = c("Core", "Accessory")))%>% mutate(Sample = "Brancorsini")

coreHet <-vcfCorsini %>% filter(Status == "Core") %>% pull(HeteroFrac)
accessHet <-vcfCorsini %>% filter(Status == "Accessory") %>% pull(HeteroFrac)
t.test(coreHet, accessHet)

# Now for Geridu
vcf <- VCFParsing("HetData/KayHet.vcf") %>% select(CHROM,FILTER, QUAL, GT)

vcfHet <- vcf %>% mutate(Hetero = ifelse(grepl("0/1|1/0", GT), T, F))

geriduGenes <- geneLengths %>% filter(Gene %in% Geridu$Gene)
vcfGeridu <- vcfHet %>% filter(!is.na(FILTER),QUAL >= 30) %>%
       	group_by(CHROM) %>% summarize(Hetero = sum(Hetero)) %>% right_join(corsiniGenes, by = c("CHROM" = "Gene")) %>%
	mutate(Hetero = replace(Hetero, is.na(Hetero),0)) %>% group_by(CHROM) %>% summarize(HeteroFrac = Hetero/Length) %>%
	mutate(HeteroFrac = replace(HeteroFrac, is.infinite(HeteroFrac), 0), Status = ifelse(CHROM %in% coreGenes, "Core", "Accessory")) %>%
	mutate(Status = factor(Status, levels = c("Core", "Accessory"))) %>% mutate(Sample = "Geridu")

coreHet <-vcfGeridu %>% filter(Status == "Core") %>% pull(HeteroFrac)
accessHet <-vcfGeridu %>% filter(Status == "Accessory") %>% pull(HeteroFrac)
t.test(coreHet, accessHet)

# Summary Statistics of both
vcfPlot <- vcfCorsini %>% bind_rows(vcfGeridu)
vcfPlot %>% group_by(Sample,Status) %>%
       	summarize(Mean = mean(HeteroFrac), SD = sd(HeteroFrac), confInt = qnorm(0.975)*SD/sqrt(length(HeteroFrac)), 
			 high = Mean + confInt, low = Mean - confInt) %>% as.data.frame()

coreHetNo <- vcfPlot %>% filter(Sample == "Brancorsini") %>% filter(HeteroFrac == 0, Status == "Core") %>%
	summarize(Genes =length(CHROM))
accessHetNo <- vcfPlot  %>% filter(Sample == "Brancorsini") %>% filter(HeteroFrac == 0, Status == "Accessory") %>%
	summarize(Genes =length(CHROM))

hetHist <- vcfPlot %>% filter(Sample == "Brancorsini") %>% ggplot(aes(x = HeteroFrac, fill = Status)) +
       	geom_histogram(position = "identity", alpha = 0.75, colour = "black") +
       #	geom_vline(xintercept = vcfPlot %>% filter(HeteroFrac < 0) %>% summarize(mean(HeteroFrac)) %>% pull(), colour = "red", lty = 2) +
       	theme_bw() +
	scale_x_log10() + annotation_logticks(sides = "b") + scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
	ylab("Genes") + xlab("P(Heterozygous)") +
	scale_fill_manual(values = c(Core = "#f8333c", Accessory = "#007dba")) +
	theme(legend.position = "bottom") +
	annotate(geom = "text", label = paste0("Core Genes with No Heterozygous SNPs: ",coreHetNo$Genes), x = 3e-3, y = 14, colour = "#f8333c") +
	annotate(geom = "text", label = paste0("Accessory Genes with No Heterozygous SNPs: ",accessHetNo$Genes), x = 3e-3, y = 12, colour = "#007dba")

# Now let's get the copy numbers involved
ancientOnly <- ancientOnly %>% left_join(vcfPlot %>% select(-Status), by = c("Sample", "Gene" = "CHROM"))
       	
model <- lm(data = ancientOnly %>% filter(HeteroFrac > 0, Sample == "Brancorsini"), log10(HeteroFrac) ~ CopyNumber*Status)
r2 <- round(summary(model)$adj.r.squared,2)
tmp <- summary(model)$fstatistic
pval <- round(pf(tmp[1],tmp[2],tmp[3], lower.tail = F),3)

copyHet <- ancientOnly %>% filter(Sample == "Brancorsini") %>% ggplot(aes(x = CopyNumber, y = HeteroFrac, colour = Status)) +
	geom_point() +
	theme_bw() +
	geom_smooth(method = "lm") +
	scale_y_log10() + annotation_logticks(sides = "l") + scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
	xlab("Copy Number") + ylab("P(Heterozygous)") +
	scale_colour_manual(values = c(Core = "#f8333c", Accessory = "#007dba")) +
	annotate(geom = "text", x = 1.6, y = 3e-3, label = bquote(R[adj]^2 == .(r2))) +
	annotate(geom = "text", x = 1.6, y = 2e-3, label = bquote(P %~~% .(pval))) +
	theme(legend.position = "bottom") 

ggarrange(hetHist, copyHet, ncol = 1, common.legend = T, legend = "bottom", align = "hv", labels = "AUTO")
ggsave(file = "Heterozygosity.pdf", height = 12, width = 9)

