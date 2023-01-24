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

ann_colors <- list(ST = c("5" = colour[3], "7" = colour[15], "8" = colour[10], "9" = colour[9], "10" = colour[6], "11" = colour[17], "12" = colour[8], "39" = colour[13], "40" = colour[2], "41" = colour[16], "42" = colour[7], "43" = colour[4], "71" = colour[11], "88" = colour[12], "102" = colour[19], "Geridu" = colour[5],"Brancorsini" = colour[1], "NF" = colour[20], "NIPH" = colour[22]))
clusterColours <- c("Western Mediterranean" = "#4D9DE0", "Fertile Crescent" = "#E15554", "Africa/America" = "#E1BC29", "Indo-Pacific" = "#3BB273", "Russia" = "#7768AE")

# Let's colour these based on their PhyloGroups 
STData <- read.delim("MLSTResults.txt")[,1:2]
rownames(STData) <- STData$Sample
STData["Brancorsini",] <- list(Sample = "Brancorsini", ST = "Brancorsini")
STData["Reference",] <- list(Sample = "Reference", ST = "7")
STData["Nodule1_S1",] <- list(Sample = "Nodule1_S1", ST = "Geridu")

STData <- STData %>% mutate(ST = gsub("[^[:alnum:] ]","", ST)) %>%  mutate(ST = replace(ST, grepl("^NIPH-*|NI_2007", Sample) & grepl("8",ST), "NIPH")) # This will highlight those odd ST 8 genomes

# Now to add the Cluster information
clusterInfo <- read.delim("../PanGenomeAnalysisFixedMapping/GenomeAccessClustered.tab")
clusterInfo[which(clusterInfo$Genome == "Geridu"),"Genome"] <- "Nodule1_S1"
rownames(clusterInfo) <- clusterInfo$Genome

# Now let's start with the Global Plot
datGlob <- as_tibble(read.delim("TempestTables/GlobalTempest.tab")) %>% filter(date > 0) %>% mutate(Phylogeny = "Global") %>% left_join(STData, by = c("tip" = "Sample"))

ann_colors$ST <- ann_colors$ST[names(ann_colors$ST) %in% datGlob$ST]

tempest <- datGlob %>%  ggplot(aes(x = date, y  = distance, group = Phylogeny)) +
       	geom_smooth(method = "lm", show.legend = F, colour = "black") +
	geom_point(alpha = 0.75, aes(colour = ST), show.legend = T) + theme_classic() + ylab("Root to Tip Divergence") + xlab("Year") +
	scale_color_manual(values = ann_colors$ST, name = "Sequence Type") + 
	theme(legend.position = "bottom") +
	guides(colour = guide_legend(nrow = 2, title.position = "top", title.hjust = 0.5))

modelGlobal <- lm(distance ~date, data = datGlob) 
r2 <- round(summary(modelGlobal)$adj.r.squared,3)
tmp <- summary(modelGlobal)$fstatistic
pval <- round(pf(tmp[1],tmp[2],tmp[3], lower.tail = F),3)

tempest <- tempest + annotate(geom = "text", x = 1700, y = 0.007,label = bquote(R[adj]^2 == .(r2))) +
	scale_x_continuous(breaks = scales:::pretty_breaks(n = 10)) +
	scale_y_continuous(breaks = scales:::pretty_breaks(n = 10)) +
	annotate(geom = "text", x = 1700, y = 0.0068, label = bquote(P == .(pval))) 

ggsave(tempest, file = "GlobalTempest.pdf", width = 6, height = 4)

# Now to do all the separate ones
files <- list.files("TempestTables", full.names = T)
files <- files[!grepl("Global", files)]

sepDat <- lapply(files ,function(x){
	location <- gsub("TempestTables\\/|Tempest.tab", "", x)
		as_tibble(read.delim(x)) %>% filter(date > 0) %>% mutate(Phylogeny = location) %>% left_join(STData, by = c("tip" = "Sample")) %>%
			left_join(clusterInfo, by = c("tip" = "Genome")) %>%
			mutate(clusters = ifelse(clusters == 1, "Western Mediterranean",
						 ifelse(clusters == 2, "Fertile Crescent",
						       	ifelse(clusters == 3, "Africa/America", ifelse(clusters == 4, "Indo-Pacific", "Russia")))))
	})

tempPlots <- lapply(sepDat, function(x){
	# Getting the location of the annotations ready
	xloc <- diff(range(x$date))/2 + range(x$date)[1]
	yloc <- diff(range(x$distance))/2 + range(x$distance)[1]

	# Making the plot
	tempPlot <- x %>%  ggplot(aes(x = date, y  = distance, group = Phylogeny)) +
	       	geom_smooth(method = "lm", show.legend = F, colour = "black") +
		geom_point(alpha = 0.75, aes(colour = ST), show.legend = T) + theme_classic() + ylab("Root to Tip Divergence") + xlab("Year") +
		scale_color_manual(values = ann_colors$ST, name = "Sequence Type") + 
		theme(legend.position = "bottom") +
		guides(colour = guide_legend(nrow = 2, title.position = "top", title.hjust = 0.5))
	# Now to calculated the model

	model <- lm(distance ~date, data = x) 
	r2 <- round(summary(model)$adj.r.squared,3)
	tmp <- summary(model)$fstatistic
	pval <- round(pf(tmp[1],tmp[2],tmp[3], lower.tail = F),3)

	tempPlot <- tempPlot + annotate(geom = "text", x = xloc, y = yloc, label = bquote(R[adj]^2 == .(r2))) +
		annotate(geom = "text", x = xloc, y = yloc * 0.99, label = bquote(P == .(pval))) 
		#scale_x_continuous(breaks = scales:::pretty_breaks(n = 10)) +
		#scale_y_continuous(breaks = scales:::pretty_breaks(n = 10)) 

})

ggarrange(plotlist = tempPlots, nrow = 2, ncol = 2, common.legend = T, legend = "bottom", labels = "auto")
ggsave("IndividualTempestPlots.pdf", width = 9, height = 6)

#################################
### Now to label via clusters ###
#################################
datGlob <- datGlob %>% left_join(clusterInfo, by = c("tip" = "Genome")) %>%
	mutate(clusters = ifelse(clusters == 1, "Western Mediterranean",
				 ifelse(clusters == 2, "Fertile Crescent",
				       	ifelse(clusters == 3, "Africa/America", ifelse(clusters == 4, "Indo-Pacific", "Russia")))))


tempest <- datGlob %>%  ggplot(aes(x = date, y  = distance, group = Phylogeny)) +
       	geom_smooth(method = "lm", show.legend = F, colour = "black") +
	geom_point(alpha = 0.75, aes(colour = clusters), show.legend = T) + theme_classic() + ylab("Root to Tip Divergence") + xlab("Year") +
	scale_color_manual(values = clusterColours, name = "Accessory PCoA Clustering") + 
	theme(legend.position = "bottom") +
	guides(colour = guide_legend(nrow = 2, title.position = "top", title.hjust = 0.5))

modelGlobal <- lm(distance ~date, data = datGlob) 
r2 <- round(summary(modelGlobal)$adj.r.squared,3)
tmp <- summary(modelGlobal)$fstatistic
pval <- round(pf(tmp[1],tmp[2],tmp[3], lower.tail = F),3)

tempest <- tempest + annotate(geom = "text", x = 1700, y = 0.007,label = bquote(R[adj]^2 == .(r2))) +
	scale_x_continuous(breaks = scales:::pretty_breaks(n = 10)) +
	scale_y_continuous(breaks = scales:::pretty_breaks(n = 10)) +
	annotate(geom = "text", x = 1700, y = 0.0068, label = bquote(P == .(pval))) 
tempest

# The separate plots
tempPlots <- lapply(sepDat, function(x){
	# Getting the location of the annotations ready
	xloc <- diff(range(x$date))/2 + range(x$date)[1]
	yloc <- diff(range(x$distance))/2 + range(x$distance)[1]

	# Making the plot
	tempPlot <- x %>% ggplot(aes(x = date, y  = distance, group = Phylogeny)) +
	       	geom_smooth(method = "lm", show.legend = F, colour = "black") +
		geom_point(alpha = 0.75, aes(colour = clusters), show.legend = T) + theme_classic() + ylab("Root to Tip Divergence") + xlab("Year") +
		scale_color_manual(values = clusterColours, name = "Accessory PCoA Clustering") + 
		theme(legend.position = "bottom") +
		guides(colour = guide_legend(nrow = 2, title.position = "top", title.hjust = 0.5))
	# Now to calculated the model

	model <- lm(distance ~date, data = x) 
	r2 <- round(summary(model)$adj.r.squared,3)
	tmp <- summary(model)$fstatistic
	pval <- round(pf(tmp[1],tmp[2],tmp[3], lower.tail = F),3)

	tempPlot <- tempPlot + annotate(geom = "text", x = xloc, y = yloc, label = bquote(R[adj]^2 == .(r2))) +
		annotate(geom = "text", x = xloc, y = yloc * 0.99, label = bquote(P == .(pval))) 
		#scale_x_continuous(breaks = scales:::pretty_breaks(n = 10)) +
		#scale_y_continuous(breaks = scales:::pretty_breaks(n = 10)) 

})
ggarrange(plotlist = tempPlots, nrow = 2, ncol = 2, common.legend = T, legend = "bottom", labels = "auto")
#ggsave("IndividualTempestPlots.pdf", width = 9, height = 6)
