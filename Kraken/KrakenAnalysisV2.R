######################################################################
# Getting the functions and Loading the libraries
library(dplyr)
library(tidyr)
library(purrr)
library(reshape2)
library(ggplot2)
library(ggpubr)
#library(ggforce)
library(ggrepel)

colour <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6',
'#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3',
'#808000', '#ffd8b1', '#000075', '#808080', '#f0f0f0', '#000000')

KrakenParsingProp <- function(fileName, Rank){
	out <- tryCatch(read.delim(fileName, header = F, comment.char = "#"), error = function(e) e)
	if(any(class(out) == "error")){
		return(data.frame("Taxon" = "Filtered", "Reads" = 0))
	}
	colnames(out) <- c("Percentage", "Reads Assigned to Clade", "Reads Assigned Directly", "Taxon Unit", "NCBI Taxon", "Taxon")

	# Filtering for requested level ID only
	profile_Species <- out[grepl(x = out$`Taxon Unit`,Rank),]

	# Fixing the names
	profile_Species$Taxon <- gsub(pattern = "\\s{2,}|-", replacement =  "", x = profile_Species$Taxon)

	# Creating the overall dataframe
	return(profile_Species[,c(6,1)])
}

KrakenParsingCount <- function(fileName, Rank){
	out <- tryCatch(read.delim(fileName, header = F, comment.char = "#"), error = function(e) e)
	if(any(class(out) == "error")){
		return(data.frame("Taxon" = "Filtered", "Reads" = 0))
	}
	colnames(out) <- c("Percentage", "Reads Assigned to Clade", "Reads Assigned Directly", "Taxon Unit", "NCBI Taxon", "Taxon")

	# Filtering for requested level ID only
	profile_Species <- out[grepl(x = out$`Taxon Unit`,Rank),]

	# Fixing the names
	profile_Species$Taxon <- gsub(pattern = "\\s{2,}|-", replacement =  "", x = profile_Species$Taxon)

	# Creating the overall dataframe
	return(profile_Species[,c(6,2)])
}

KrakenAnalysis <- function(profile_grouped, noise = 0.01){
	# If NA, set to 0
	profile_grouped[,-1] <- as.data.frame(sapply(profile_grouped[,-1], function(x){ifelse(is.na(x), 0, x)}))
	profileZeroes <-  profile_grouped
	
	# Now to remove those taxa which account for 1% of the reads in each column
	unknown = list("Other")
	for(i in 2:ncol(profile_grouped)){
		index <- profile_grouped[,i] < noise
		unknown[i] <- sum(profile_grouped[index,i], na.rm = T)
		profile_grouped[index,i] <- 0
	}

	profile_grouped <- profile_grouped[rowSums(profile_grouped[,-1]) > 0,]
	profile_grouped[nrow(profile_grouped) + 1,] <- unknown
	digestFreq <- profile_grouped %>% group_by(Taxon) %>% pivot_longer(c(everything(), -Taxon), names_to = "Sample", values_to = "Reads") %>%
		ungroup()
#	finalRaw <- profile_grouped %>% group_by(Taxon) %>% pivot_longer(c(everything(), -Taxon), names_to = "Sample", values_to = "Reads") %>%
#		ungroup()
	
	return(digestFreq)
	#return(list("Proportional" = digestFreq, "Reads" = finalRaw, "OrigTable" = profileZeroes))
	
}

theme_set(theme_classic())
############
###Kraken###
############
# Getting the files of interest
report_files <- list.files(path ="CombinedReports", full.names = T)
ranked = "^F$"

# Getting the ID numbers - Used to ease IDing the files of interest
#num <- gsub(".tab","", basename(report_files))
num <- c("Brancorsini - Nodule", "Brancorsini - Extraction Blank", "Geridu")

tmp <- lapply(report_files, function(x){KrakenParsingProp(x,ranked)})
krakentable <- as_tibble(reduce(tmp, by = "Taxon",full_join))
tmp <- lapply(report_files, function(x){KrakenParsingCount(x,ranked)})
krakentableCount <- as_tibble(reduce(tmp, by = "Taxon",full_join))
colnames(krakentable)[-1] <- colnames(krakentableCount)[-1] <-  num

krakentable <- krakentable %>% pivot_longer(c(everything(), -Taxon)) %>% group_by(name) %>% mutate(value = value/sum(value, na.rm = T)) %>% 
	pivot_wider(Taxon)
krakenProp <- KrakenAnalysis(krakentable, 0.01)

########################################
### Preparing the Proportional Data ####
########################################
#Sample Translations
krakenabund <- krakentableCount %>% filter(Taxon %in% unique(krakenProp$Taxon)) %>% pivot_longer(c(everything(), -Taxon), names_to = "Sample") %>% group_by(Sample) %>%
	summarize(Abundance = sum(value, na.rm = T))


# What are we looking at specifically
tmp <- krakenProp %>% pull(Taxon) %>% unique() %>% sort()
ind <- c(which(tmp == "Hominidae"), which(tmp == "Other"))
ind2 <- which(tmp == "Brucellaceae")
ind3 <- c(ind, ind2)
#plotDf <- krakenProp %>%
#       	mutate(Taxon = factor(Taxon, levels = c(tmp[ind2],tmp[-ind3],tmp[ind])), Sample = factor(Sample, levels = c("Brancorsini", "Geridu", "Library Blank", "Extraction Blank")))# %>% mutate(Type = unlist(GroupList[Sample]))
plotDf <- krakenProp %>%
       	mutate(Taxon = factor(Taxon, levels = c(tmp[ind2],tmp[-ind3],tmp[ind])))#, Sample = factor(Sample, levels = c("Brancorsini", "Geridu", "Water Blank")))# %>% mutate(Type = unlist(GroupList[Sample]))

fam <- plotDf %>% 
	ggplot(aes(x = Sample, y = Reads, fill = Taxon)) +
	geom_col() + scale_fill_manual(values = colour) +
	theme(legend.text = element_text(face = "italic")) +
	#theme(axis.text.x = element_blank(),axis.title.x = element_blank(),legend.text = element_text(face = "italic")) +
	scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
	geom_text(data = krakenabund, inherit.aes = F, aes(x = Sample, label = Abundance, y = 1.02), size = 2.5) +
	guides(fill = guide_legend(title = "Family", ncol = 2)) + ylab("Abundance")
ggsave(fam,file = "Kraken2Family0.012.pdf", width = 9, height = 6)
fam

#############
### Genus ###
#############
ranked = "^G$"
tmp <- lapply(report_files, function(x){KrakenParsingProp(x,ranked)})
krakentable <- as_tibble(reduce(tmp, by = "Taxon",full_join))
tmp <- lapply(report_files, function(x){KrakenParsingCount(x,ranked)})
krakentableCount <- as_tibble(reduce(tmp, by = "Taxon",full_join))
colnames(krakentable)[-1] <- colnames(krakentableCount)[-1] <-  num

krakentable <- krakentable %>% pivot_longer(c(everything(), -Taxon)) %>% group_by(name) %>% mutate(value = value/sum(value, na.rm = T)) %>% 
	pivot_wider(Taxon)
krakenProp <- KrakenAnalysis(krakentable, 0.01)

########################################
### Preparing the Proportional Data ####
########################################
krakenabund <- krakentableCount %>% filter(Taxon %in% unique(krakenProp$Taxon)) %>% pivot_longer(c(everything(), -Taxon), names_to = "Sample") %>% group_by(Sample) %>%
	summarize(Abundance = sum(value, na.rm = T))

# What are we looking at specifically
tmp <- krakenProp %>% pull(Taxon) %>% unique() %>% sort()
ind <- c(which(tmp == "Homo"), which(tmp == "Other"))
ind2 <- which(tmp == "Brucella")
ind3 <- c(ind, ind2)
plotDf <- krakenProp %>%
       	mutate(Taxon = factor(Taxon, levels = c(tmp[ind2],tmp[-ind3],tmp[ind]))) %>%
	mutate(Sample = factor(Sample, levels = c("Brancorsini - Nodule", "Brancorsini - Extraction Blank", "Geridu")))#, Sample = factor(Sample, levels = c("Brancorsini", "Geridu", "Water Blank")))# %>% mutate(Type = unlist(GroupList[Sample]))
#plotDf <- krakenProp %>%
#       	mutate(Taxon = factor(Taxon, levels = c(tmp[ind2],tmp[-ind3],tmp[ind])), Sample = factor(Sample, levels = c("Brancorsini", "Geridu", "Library Blank", "Extraction Blank")))# %>% mutate(Type = unlist(GroupList[Sample]))
#
gen <- plotDf %>%
	ggplot(aes(x = Sample, y = Reads, fill = Taxon)) +
	geom_col() + scale_fill_manual(values = colour) +
	theme(legend.text = element_text(face = "italic")) +
	#theme(axis.text.x = element_blank(),axis.title.x = element_blank(),legend.text = element_text(face = "italic")) +
	scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
	geom_text(data = krakenabund, inherit.aes = F, aes(x = Sample, label = Abundance, y = 1.02), size = 2.5) +
	guides(fill = guide_legend(title = "Genus", ncol = 4, title.position = "top", title.hjust = 0.5)) +
       	ylab("Abundance") + theme(legend.position = "bottom")
gen
ggsave(gen,file = "Kraken2Genus.pdf", width = 6, height = 4)
###############
### Species ###
###############
ranked = "^S$"
tmp <- lapply(report_files, function(x){KrakenParsingProp(x,ranked)})
krakentable <- as_tibble(reduce(tmp, by = "Taxon",full_join))
tmp <- lapply(report_files, function(x){KrakenParsingCount(x,ranked)})
krakentableCount <- as_tibble(reduce(tmp, by = "Taxon",full_join))
colnames(krakentable)[-1] <- colnames(krakentableCount)[-1] <-  num

krakentable <- krakentable %>% pivot_longer(c(everything(), -Taxon)) %>% group_by(name) %>% mutate(value = value/sum(value, na.rm = T)) %>% 
	pivot_wider(Taxon)
krakenProp <- KrakenAnalysis(krakentable, 0.01)

########################################
### Preparing the Proportional Data ####
########################################
krakenabund <- krakentableCount %>% filter(Taxon %in% unique(krakenProp$Taxon)) %>% pivot_longer(c(everything(), -Taxon), names_to = "Sample") %>% group_by(Sample) %>%
	summarize(Abundance = sum(value, na.rm = T))
# What are we looking at specifically
tmp <- krakenProp %>% pull(Taxon) %>% unique() %>% sort()
ind <- c(which(tmp == "Homo sapiens"), which(tmp == "Other"))
ind2 <- which(tmp == "Brucella melitensis")
ind3 <- c(ind, ind2)
plotDf <- krakenProp %>%
       	mutate(Taxon = factor(Taxon, levels = c(tmp[ind2],tmp[-ind3],tmp[ind])), Sample = factor(Sample, levels = c("Brancorsini", "Geridu", "Water Blank")))# %>% mutate(Type = unlist(GroupList[Sample]))
#plotDf <- krakenProp %>%
#       	mutate(Taxon = factor(Taxon, levels = c(tmp[ind2],tmp[-ind3],tmp[ind])), Sample = factor(Sample, levels = c("Brancorsini", "Geridu", "Library Blank", "Extraction Blank")))# %>% mutate(Type = unlist(GroupList[Sample]))
#
spe <- plotDf %>%
	ggplot(aes(x = Sample, y = Reads, fill = Taxon)) +
	geom_col() + scale_fill_manual(values = colour) +
	theme(legend.text = element_text(face = "italic")) +
	#theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),legend.text = element_text(face = "italic")) +
	scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
	geom_text(data = krakenabund, inherit.aes = F, aes(x = Sample, label = Abundance, y = 1.02), size = 2.5) +
	guides(fill = guide_legend(title = "Species", ncol = 2)) + ylab("Abundance")
spe
ggsave("KrakenSpecies01.pdf", width = 9, height = 6)

ggarrange(fam, gen,spe, ncol = 1, align = "hv", labels = "AUTO")
ggsave("KrakenNoUnclass.png", height = 12, width = 9)
#spe
#ggsave("../SpeciesKrakenV2.pdf", width = 8, height = 6)
