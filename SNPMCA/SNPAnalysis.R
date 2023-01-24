library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(ggExtra)
library(scales)
library(cluster)
library(pbapply)
library(parallel)
library(purrr)
library(reshape2)
library(FactoMineR)
library(gtools)
library(pheatmap)
library(factoextra)
library(emmeans)
#library(ggvenn)
#library(pvclust)
#library(ggdendro)
theme_set(theme_classic())

## Functions

SNPMCAPlotting <- function(snpsMatrix, metadata, baryVar = NULL){
	mcaResults <- MCA(snpsMatrix, graph = F)
	eigenVals <- get_eigenvalue(mcaResults) %>% as_tibble()
	
	mcaCoord <- mcaResults$ind$coord %>% as.data.frame()
	colnames(mcaCoord) <- gsub(" ","", colnames(mcaCoord))
	mcaCoord$Sample <- rownames(mcaCoord)
	mcaCoord <- as_tibble(mcaCoord) %>% left_join(metadata) %>% 
		mutate(ST = as.character(ST),
		       ST = replace(ST, Sample == "Brancorsini", "Brancorsini"),
		       ST = replace(ST, Sample == "Nodule1_S1", "Geridu"))
	
	# Calculating the Barycenters
	if(!is.null(baryVar)){
		if(baryVar %in% colnames(mcaCoord)){
			barycenter <- mcaCoord %>% group_by(.data[[baryVar]]) %>% select(starts_with("Dim")) %>%
				summarise_all(mean)
			plotBary <- T
		}else{
			stop("Grouping Variable for Barycentre not found!")
		}
	}else{
		baryVar = "Host"
		plotBary <- F
	}
	#fviz_nbclust(x = mcaCoord[,1:5], method = "silhouette", FUNcluster = clara)
	#mcaCoord$Clusters <- as.factor(clara(mcaCoord[,1:5], k = 7)$clustering)
		
	mcaInd <- mcaCoord %>% ggplot(aes(x = Dim1, y = Dim2, fill = .data[[baryVar]], colour = .data[[baryVar]])) +
		geom_vline(xintercept = 0, colour = "grey90", lty = 2) +
		geom_hline(yintercept = 0, colour = "grey90", lty = 2) +
		stat_ellipse(level = 0.95, alpha = 0.5) +
		geom_point(alpha = ifelse(plotBary, 0.5,1)) +
		xlab(paste0("Dim 1 (",round(eigenVals$variance.percent[1], 3),"%)")) +
		ylab(paste0("Dim 2 (",round(eigenVals$variance.percent[2], 3),"%)")) +
		scale_colour_manual(values = annColours[[baryVar]]) +
		scale_fill_manual(values = annColours[[baryVar]]) +
		guides(colour = guide_legend(nrow = 2), fill = guide_legend(nrow = 2))
	if(plotBary) mcaInd <- mcaInd + geom_point(size = 2, shape = "square", data = barycenter)
	
	mcaIndpt2 <- mcaCoord %>% ggplot(aes(x = Dim1, y = Dim3, colour = .data[[baryVar]])) +
		geom_vline(xintercept = 0, colour = "grey90", lty = 2) +
		geom_hline(yintercept = 0, colour = "grey90", lty = 2) +
		stat_ellipse(level = 0.95, alpha = 0.5) +
		geom_point(alpha = ifelse(plotBary, 0.5,1)) +
		xlab(paste0("Dim 1 (",round(eigenVals$variance.percent[1], 3),"%)")) +
		ylab(paste0("Dim 3 (",round(eigenVals$variance.percent[3], 3),"%)")) +
		scale_colour_manual(values = annColours[[baryVar]]) +
		scale_fill_manual(values = annColours[[baryVar]]) 
	if(plotBary) mcaIndpt2 <- mcaIndpt2 + geom_point(size = 2, shape = "square", data = barycenter)
	
	# Now looking at the variables
	mcaCoordVar <- mcaResults$var$coord %>% as.data.frame()
	colnames(mcaCoordVar) <- gsub(" ","", colnames(mcaCoordVar))
	mcaCoordVar$SNP <- rownames(mcaCoordVar)
	
	mcaCoordVarCos <- mcaResults$var$cos2 %>% as.data.frame()
	colnames(mcaCoordVarCos) <- gsub(" ","", colnames(mcaCoordVarCos))
	mcaCoordCosDim1Dim2 <- mcaCoordVarCos %>% select(Dim1, Dim2) %>% rowSums()
	
	mcaCoordVar <- mcaCoordVar %>% left_join(data.frame(SNP = names(mcaCoordCosDim1Dim2), Cos2 = mcaCoordCosDim1Dim2))
	mcaCoordVar <- mcaCoordVar %>% mutate(Cos2Imp = ifelse(Cos2 >= 0.1, "Important", "Unimportant")) %>%
		mutate(Label = SNP, Label = replace(Label, Cos2 < 0.1, NA)) %>%
			mutate(SNPType = ifelse(grepl("_X$", SNP), "Gene Absent",
						ifelse(grepl("_N$", SNP), "Ambiguous Base",
						       ifelse(grepl("_S$", SNP), "Reference",
							      ifelse(grepl("_M$", SNP), "No Coverage","SNP")))))
	
	mcaVarpt1 <- mcaCoordVar %>% ggplot(aes(x = Dim1, y = Dim2, colour = SNPType)) +
		#geom_rect(aes(xmin = -0.5, xmax = 0.5, ymin = -0.5, ymax = 0.5), fill = NA, colour = "black", lty = 2) +
		geom_vline(xintercept = 0, colour = "grey90", lty = 2) +
		geom_hline(yintercept = 0, colour = "grey90", lty = 2) +
		#geom_point(colour = ifelse(mcaCoordVar$Cos2 >= 0.1, "red", "black")) +
		geom_point(alpha = 0.5) +
		#scale_colour_manual(values = c("Important" = colour[1], "Unimportant" = "grey50"), name = "Cos2 Importance") +
		scale_colour_manual(values = annColours$SNPType, name = "SNP") +
		#scale_colour_viridis_c() +
		xlab(paste0("Dim 1 (",round(eigenVals$variance.percent[1], 3),"%)")) +
		ylab(paste0("Dim 2 (",round(eigenVals$variance.percent[2], 3),"%)")) 
		#geom_text_repel(show.legend = F) #+ theme(legend.position = "bottom")
	
	mcaVarpt2 <- mcaCoordVar %>% ggplot(aes(x = Dim1, y = Dim3, colour = SNPType, label = Label)) +
		#geom_rect(aes(xmin = -0.5, xmax = 0.5, ymin = -0.5, ymax = 0.5), fill = NA, colour = "black", lty = 2) +
		geom_vline(xintercept = 0, colour = "grey90", lty = 2) +
		geom_hline(yintercept = 0, colour = "grey90", lty = 2) +
		#geom_point(colour = ifelse(mcaCoordVar$Cos2 >= 0.1, "red", "black")) +
		geom_point(alpha = 0.5) +
		#geom_density_2d() +
		#scale_colour_manual(values = c("Important" = colour[1], "Unimportant" = "grey50"), name = "Cos2 Importance") +
		scale_colour_manual(values = annColours$SNPType, name = "SNP") +
		#scale_colour_viridis_c() +
		xlab(paste0("Dim 1 (",round(eigenVals$variance.percent[1], 3),"%)")) +
		ylab(paste0("Dim 3 (",round(eigenVals$variance.percent[3], 3),"%)")) 
		#geom_text_repel(show.legend = F) #+ theme(legend.position = "bottom")
	
	# Making the plots
	varPlot <- ggarrange(mcaVarpt1, mcaVarpt2, nrow =1, align = "hv", common.legend =T, legend = "bottom", labels = c("c", "d"))
	IndPlot <- ggarrange(mcaInd, mcaIndpt2, nrow =1, align = "hv", common.legend =T, legend = "bottom", labels = "auto")
	mcaPlot <- ggarrange(IndPlot, varPlot, align = "hv", nrow = 2)
	return(list("Plots" = mcaPlot, "Data" = mcaCoord, "bary" = barycenter))
}
SNPMCAIndOnly <- function(snpsMatrix, metadata, baryVar = NULL){
	mcaResults <- MCA(snpsMatrix, graph = F)
	eigenVals <- get_eigenvalue(mcaResults) %>% as_tibble()
	
	mcaCoord <- mcaResults$ind$coord %>% as.data.frame()
	colnames(mcaCoord) <- gsub(" ","", colnames(mcaCoord))
	mcaCoord$Sample <- rownames(mcaCoord)
	mcaCoord <- as_tibble(mcaCoord) %>% left_join(metadata) %>% 
		mutate(ST = as.character(ST),
		       ST = replace(ST, Sample == "Brancorsini", "Brancorsini"),
		       ST = replace(ST, Sample == "Nodule1_S1", "Geridu"))
	
	# Calculating the Barycenters
	if(!is.null(baryVar)){
		if(baryVar %in% colnames(mcaCoord)){
			barycenter <- mcaCoord %>% group_by(.data[[baryVar]]) %>% select(starts_with("Dim")) %>%
				summarise_all(mean)
			plotBary <- T
		}else{
			stop("Grouping Variable for Barycentre not found!")
		}
	}else{
		baryVar = "Host"
		plotBary <- F
	}
	#fviz_nbclust(x = mcaCoord[,1:5], method = "silhouette", FUNcluster = clara)
	#mcaCoord$Clusters <- as.factor(clara(mcaCoord[,1:5], k = 7)$clustering)
		
	mcaInd <- mcaCoord %>% ggplot(aes(x = Dim1, y = Dim2, fill = .data[[baryVar]], colour = .data[[baryVar]])) +
		geom_vline(xintercept = 0, colour = "grey90", lty = 2) +
		geom_hline(yintercept = 0, colour = "grey90", lty = 2) +
		#stat_ellipse(level = 0.95, alpha = 0.5) +
		geom_point(alpha = ifelse(plotBary, 0.5,1)) +
		xlab(paste0("Dim 1 (",round(eigenVals$variance.percent[1], 3),"%)")) +
		ylab(paste0("Dim 2 (",round(eigenVals$variance.percent[2], 3),"%)")) +
		scale_colour_manual(values = annColours[[baryVar]]) +
		scale_fill_manual(values = annColours[[baryVar]]) +
		guides(colour = guide_legend(nrow = 2), fill = guide_legend(nrow = 2))
	#if(plotBary) mcaInd <- mcaInd + geom_point(size = 2, shape = "square", data = barycenter)
	
	mcaIndpt2 <- mcaCoord %>% ggplot(aes(x = Dim1, y = Dim3, colour = .data[[baryVar]])) +
		geom_vline(xintercept = 0, colour = "grey90", lty = 2) +
		geom_hline(yintercept = 0, colour = "grey90", lty = 2) +
		#stat_ellipse(level = 0.95, alpha = 0.5) +
		geom_point(alpha = ifelse(plotBary, 0.5,1)) +
		xlab(paste0("Dim 1 (",round(eigenVals$variance.percent[1], 3),"%)")) +
		ylab(paste0("Dim 3 (",round(eigenVals$variance.percent[3], 3),"%)")) +
		scale_colour_manual(values = annColours[[baryVar]]) +
		scale_fill_manual(values = annColours[[baryVar]]) 
	#if(plotBary) mcaIndpt2 <- mcaIndpt2 + geom_point(size = 2, shape = "square", data = barycenter)
	
	IndPlot <- ggarrange(mcaInd, mcaIndpt2, nrow =1, align = "hv", common.legend =T, legend = "bottom", labels = "auto")
	return(list("Plots" = IndPlot, "Data" = mcaCoord, "bary" = barycenter))
}
############
colour <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3',
'#808000', '#ffd8b1', '#000075', '#808080', '#f0f0f0', '#000000')

annColours <- list(ST = c("5" = colour[3], "7" = colour[15], "8" = colour[10], "9" = colour[9], "10" = colour[6], "11" = colour[17], "12" = colour[8], "39" = colour[13], "40" = colour[2], "41" = colour[16], "42" = colour[7], "43" = colour[4], "71" = colour[11], "88" = colour[12], "102" = colour[19],"Brancorsini" = colour[1], "NF" = colour[20], "NIPH" = colour[22]),
	Norway = c("Normal" = "#00205B" , "Odd" = "#BA0C2F", "Other" = "#FFFFFF"),
	SNPType = c("Reference" = "#3BB273", "SNP" = "#4D9DE0", "Ambiguous Base" = "#7768AE", "No Coverage" = "#E1BC29", "Gene Absent" = "#E15554"),
	clusters = c("Western Mediterranean" = "#4D9DE0", "Fertile Crescent" = "#E15554", "Africa/America" = "#E1BC29", "Indo-Pacific" = "#3BB273", "Russia" = "#7768AE"))
#annColours <- list(ST = c("5" = colour[3], "7" = colour[15], "8" = colour[10], "9" = colour[9], "10" = colour[6], "11" = colour[17], "12" = colour[8], "39" = colour[13], "40" = colour[2], "41" = colour[16], "42" = colour[7], "43" = colour[4], "71" = colour[11], "88" = colour[12], "102" = colour[19], "Geridu" = colour[5],"Brancorsini" = colour[1], "NF" = colour[20], "Reference" = colour[22]),
#	Norway = c("Normal" = "#00205B" , "Odd" = "#BA0C2F", "Other" = "#FFFFFF"),
#	SNPType = c("Reference" = "#3BB273", "SNP" = "#4D9DE0", "Ambiguous Base" = "#7768AE", "No Coverage" = "#E1BC29", "Gene Absent" = "#E15554"),
#	clusters = c("Western Mediterranean" = "#4D9DE0", "Fertile Crescent" = "#E15554", "Africa/America" = "#E1BC29", "Indo-Pacific" = "#3BB273", "Russia" = "#7768AE"))

# Reading in the metadata information
MLSTResults <- read.delim("../PanGenomeAnalysisFixedMapping/MLSTResults.txt") %>% as_tibble()
MLSTResults <- MLSTResults %>% select(Sample, ST) %>% mutate(ST = gsub("\\*| ","", ST))
colnames(MLSTResults)[1] <- "Genome"

# Reading the metadata table
metaData <- read.delim("../PanGenomeAnalysisFixedMapping/MetadataAll.tab", header =T) %>% as_tibble

clusterInfo <- read.delim("../PanGenomeAnalysisFixedMapping/GenomeAccessClustered.tab") %>%
	mutate(clusters = ifelse(clusters == 1, "Western Mediterranean",
				 ifelse(clusters == 2, "Fertile Crescent",
				       	ifelse(clusters == 3, "Africa/America", ifelse(clusters == 4, "Indo-Pacific", "Russia"))))) %>% as_tibble()
#clusterInfo[which(clusterInfo$Genome == "Geridu"),"Genome"] <- "Nodule1"

metaData <- metaData %>% left_join(clusterInfo, by = c("Sample" = "Genome")) 
metaData[which(metaData$Sample == "Brancorsini"),"clusters"] <- "Western Mediterranean"
metaData[which(metaData$Sample == "Brancorsini"),"Sample"] <- "Brancorsini"

#metaData[which(metaData$Sample == "KayBMel"),"clusters"] <- "Western Mediterranean"
#metaData[which(metaData$Sample == "KayBMel"),"Sample"] <- "Nodule_S1"

metaData <- metaData %>% select(Sample, ST, clusters) %>% distinct()%>%
       	mutate(ST = replace(ST, grepl("^NIPH-*|NI_2007", Sample) & grepl("8",ST), "NIPH")) %>%
	mutate(ST = replace(ST, grepl("Brancorsini", Sample), "Brancorsini")) 

# Now to read in the data
files <- list.files("SNPCSV", full.names = T)

snpList <- pblapply(cl = detectCores() - 1, files, FUN = function(x){

		 # Let's get the WELL ID out in the open
		 sampleName <- gsub(".*/|\\.csv.*","",x)

		 snps <- read.csv(x, colClasses = c("character", "numeric", rep("character",4), rep(NA, 8)))
		 if(nrow(snps) > 0 ){
		 	snps$Sample <- sampleName
		 	return(snps)
		 }else{
			 return(snps)
		 }
	})
names(snpList) <- gsub(".*/|_.*|\\.csv.*","",files)

# Binding the rows together
cat(sum(sapply(snpList, nrow) == 0), "Samples had no SNPs Detected\n")
noSNPs <- names(snpList[sapply(snpList, nrow) == 0])
snpList <- snpList[sapply(snpList, nrow) > 0] # Need to remove the dataframes with nothing in there....

# We're only interested in SNPs for now and we need to expand and contract the dataframe to figure which SNPs are missing
#snps <- bind_rows(snpList) %>% as_tibble() %>% select(CHROM, POS, TYPE,REF,ALT,Sample) %>% filter(TYPE == "snp") %>% unite(col = Location, CHROM, POS, sep = ";")

#uniquesnps <- snps %>% select(Location, REF) %>% distinct()
#
#snps <- snps %>% pivot_wider(id_cols = Location,names_from = Sample, values_from = ALT) %>%
#	pivot_longer(-matches("Location"), values_to = "ALT", names_to = "Sample") %>% 
#	left_join(uniquesnps) %>% select(Location, REF, ALT, Sample) 
#snpList <- split(snps, f = snps$Sample)
#
# Now to figure what I'm missing here by sample 
#test <- pblapply(cl = detectCores() - 1, snpList, function(x){
#	# Now to find which positions I have missing information
#	missingInfo <- x %>% filter(is.na(ALT))
#
#	# Pulling out the Sample and loading in the information needed
#	Sample <- x$Sample[1]
#
#	depths <- as_tibble(read.delim(paste0("ModernDepthsFiltered/",Sample,".tab.gz"), header = F, col.names = c("Gene", "Pos", "Coverage"))) %>%
#	        unite(Location, Gene, Pos, sep = ";") %>% filter(Location %in% missingInfo$Location) %>%
#		mutate(SNP = ifelse(Coverage >= 10, "S", ifelse(Coverage < 10 & Coverage >= 1, "N", "M")))
#
#	missingInfoReplaced <- missingInfo %>% left_join(depths, by = "Location") %>% mutate(SNP = replace(SNP, is.na(SNP), "X")) %>%
#	       	select(Location, REF, SNP, Sample) %>% rename(ALT = SNP)
#
#	filledSNPs <- x %>% filter(!is.na(ALT)) %>% bind_rows(missingInfoReplaced)
#	
#	return(filledSNPs)
#})

###################################
### Doing the SNP MCA Analysis ### 
###################################
# The simple way
snpsUpdated <- bind_rows(snpList) %>% unite(col = Location, CHROM, POS) %>% filter(TYPE == "snp")

# The more complicated way.....
# Before making the SNP matrix, need to filter out the genes not present in our P/A analysis
#ancientDepths <- read.delim("../PanGenomeAnalysisFixedMapping/AncientDepths.tab") %>% as_tibble() %>% mutate(Genome = replace(Genome, Genome == "Geridu", "Nodule1_S1"))
#
#branPresent <- ancientDepths %>% filter(CV <= 1.5, GCCorrected >= 10, Genome == "Brancorsini") %>% pull(Gene)
#gerPresent <- ancientDepths %>% filter(CV <= 1.5, GCCorrected >= 1, Genome != "Brancorsini") %>% pull(Gene)
#
#test$Brancorsini <- test$Brancorsini %>% separate(Location, into = c("Gene", "Pos"), sep = ";") %>% mutate(ALT = replace(ALT, !(Gene %in% branPresent), "X")) %>% unite(Location, Gene, Pos, sep = ";")
#test$Nodule1_S1 <- test$Nodule1_S1 %>% separate(Location, into = c("Gene", "Pos"), sep = ";") %>% mutate(ALT = replace(ALT, !(Gene %in% gerPresent), "X")) %>% unite(Location, Gene, Pos, sep = ";")
#
#snpsUpdated <- bind_rows(test) %>% mutate(Location = gsub(";","_", Location)) %>% filter(Sample != "Nodule1_S1")
#
###
snpsMatrix <- snpsUpdated %>% dcast(Location ~ Sample, value.var = "ALT", fill = "X") 
rownames(snpsMatrix) <- snpsMatrix$Location
snpsMatrix <- snpsMatrix[,-1]
snpsMatrix <- t(snpsMatrix) %>% as.data.frame()

# Want to try removing the singletons, could be strongly influencing the results
nonSingletonVariants <- apply(snpsMatrix, MARGIN = 2, FUN = function(x){
		     tmp <- table(x) %>% as.vector()
		     if((sum(tmp == 1) & length(tmp) == 2) | length(tmp) == 1){
			     return(FALSE) # Removing the pesky singletons
		     }else{
			     return(TRUE)
		     }
		  })

cat(sum(!nonSingletonVariants), "Singleton SNPs Removed\n")
snpsMatrix <- snpsMatrix[,nonSingletonVariants] # Removes 149 Variants

# Removing the samples that are just X
#snpsMatrix <- snpsMatrix[sampleNames,]

# Continuing our Programming
snpsMatrix <- snpsMatrix %>% mutate_all(as.factor)

# Trying an MCA
snpAllIndOnly <- SNPMCAIndOnly(snpsMatrix, metaData, baryVar = "clusters")
snpAllIndOnly$Plots
ggsave(snpAllIndOnly$Plots, file = "Plots/AllSNPsIndexOnly.pdf", width = 9, height = 6)

###########################################
### Now to look only at the core genome ###
###########################################
coreGenes <- read.table("../PanGenomeAnalysisFixedMapping/AssemblyOnly/CoreGenes99.list") %>% pull()

#snpscore <- bind_rows(test) %>% mutate(CHROM = gsub(";.*","", Location), Location = gsub(";","_", Location)) %>% filter(CHROM %in% coreGenes, Sample != "Nodule1_S1") 
snpscore <- bind_rows(snpList) %>% filter(CHROM %in% coreGenes, TYPE == "snp") %>% unite(col = Location, CHROM, POS) # Simple

snpsMatrix <- snpscore %>% dcast(Location ~ Sample, value.var = "ALT", fill = "X") 
rownames(snpsMatrix) <- snpsMatrix$Location
snpsMatrix <- snpsMatrix[,-1]
snpsMatrix <- t(snpsMatrix) %>% as.data.frame()

# Want to try removing the singletons, could be strongly influencing the results
nonSingletonVariants <- apply(snpsMatrix, MARGIN = 2, FUN = function(x){
		     tmp <- table(x) %>% as.vector()
		     if((sum(tmp == 1) & length(tmp) == 2) | length(tmp) == 1){
			     return(FALSE) # Removing the pesky singletons
		     }else{
			     return(TRUE)
		     }
		  })

cat(sum(!nonSingletonVariants), "Singleton SNPs Removed\n")
snpsMatrix <- snpsMatrix[,nonSingletonVariants] # Removes 149 Variants

# Removing the samples that are just X
#snpsMatrix <- snpsMatrix[sampleNames,]

# Continuing our Programming
snpsMatrix <- snpsMatrix %>% mutate_all(as.factor)

# Trying an MCA
snpCore <- SNPMCAIndOnly(snpsMatrix, metaData, baryVar = "clusters")
snpCore$Plots
snpCoreST <- SNPMCAIndOnly(snpsMatrix, metaData, baryVar = "ST")
snpCoreST$Plots
ggsave(snpCoreST$Plots, file = "PlotsAssembly/CoreSNP.pdf", width = 9, height = 6)

################################################
### Now to look only at the Accessory genome ###
################################################
#snpsaccess <- bind_rows(test) %>% mutate(CHROM = gsub(";.*","", Location), Location = gsub(";","_", Location)) %>% filter(!(CHROM %in% coreGenes), Sample != "Nodule1_S1")
snpsaccess <- bind_rows(snpList) %>% filter(!(CHROM %in% coreGenes), TYPE == "snp") %>% unite(col = Location, CHROM, POS)

snpsMatrix <- snpsaccess %>% dcast(Location ~ Sample, value.var = "ALT", fill = "X") 
rownames(snpsMatrix) <- snpsMatrix$Location
snpsMatrix <- snpsMatrix[,-1]
snpsMatrix <- t(snpsMatrix) %>% as.data.frame()

# Want to try removing the singletons, could be strongly influencing the results
nonSingletonVariants <- apply(snpsMatrix, MARGIN = 2, FUN = function(x){
		     tmp <- table(x) %>% as.vector()
		     if((sum(tmp == 1) & length(tmp) == 2) | length(tmp) == 1){
			     return(FALSE) # Removing the pesky singletons
		     }else{
			     return(TRUE)
		     }
		  })

cat(sum(!nonSingletonVariants), "Singleton SNPs Removed\n")
snpsMatrix <- snpsMatrix[,nonSingletonVariants] # Removes 149 Variants

# Removing the samples that are just X
#snpsMatrix <- snpsMatrix[sampleNames,]

# Continuing our Programming
snpsMatrix <- snpsMatrix %>% mutate_all(as.factor)

# Trying an MCA
snpaccess <- SNPMCAIndOnly(snpsMatrix, metaData, baryVar = "clusters")
ggsave(snpaccess$Plots, file = "PlotsAssembly/AccessorySNPClustered.pdf", width = 9, height = 6)
snpaccessST <- SNPMCAIndOnly(snpsMatrix, metaData, baryVar = "ST")
ggsave(snpaccessST$Plots, file = "PlotsAssembly/AccessorySNPST.pdf", width = 9, height = 6)
##########################################################
### Finally, limiting ourselves to the Virulence Genes ###
##########################################################
virGenes <- read.table("../PanGenomeAnalysisFixedMapping/IdentifiedVirulence.list") %>% pull()

#snpsvir <- bind_rows(test) %>% mutate(CHROM = gsub(";.*","", Location), Location = gsub(";","_", Location)) %>% filter(CHROM %in% virGenes, Sample != "Nodule1_S1")
snpsvir <- bind_rows(snpList) %>% filter(CHROM %in% virGenes, TYPE == "snp") %>% unite(col = Location, CHROM, POS) # Simple way

snpsMatrix <- snpsvir %>% dcast(Location ~ Sample, value.var = "ALT", fill = "X") 
rownames(snpsMatrix) <- snpsMatrix$Location
snpsMatrix <- snpsMatrix[,-1]
snpsMatrix <- t(snpsMatrix) %>% as.data.frame()

# Want to try removing the singletons, could be strongly influencing the results
nonSingletonVariants <- apply(snpsMatrix, MARGIN = 2, FUN = function(x){
		     tmp <- table(x) %>% as.vector()
		     if((sum(tmp == 1) & length(tmp) == 2) | length(tmp) == 1){
			     return(FALSE) # Removing the pesky singletons
		     }else{
			     return(TRUE)
		     }
		  })

cat(sum(!nonSingletonVariants), "Singleton SNPs Removed\n")
snpsMatrix <- snpsMatrix[,nonSingletonVariants] # Removes 149 Variants

# Removing the samples that are just X
#snpsMatrix <- snpsMatrix[sampleNames,]

# Continuing our Programming
snpsMatrix <- snpsMatrix %>% mutate_all(as.factor)

# Trying an MCA
snpvir <- SNPMCAIndOnly(snpsMatrix, metaData, baryVar = "clusters")
snpvir$Plots
ggsave(snpvir$Plots, file = "Plots/VirulenceSNPsIndexOnly.pdf", width = 9, height = 6)
