# Loading Libraries, functions, and variables that will be needed
library(dplyr)
library(tidyr)
library(scales)
library(purrr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(cowplot)
library(MASS)
library(lme4)
library(emmeans)
library(gamlss)

RateParsing <- function(fileName){
	tallyName <- gsub(".*\\/","",gsub("\\/Stats_out_.*","", fileName))
	tallyName <- gsub("_.*", "", tallyName)
	test <- read.csv(fileName)
	tmp <- tibble("Sample" = tallyName, "DeltaS" = test$DeltaS[1], "Std" = test$DeltaS[2])
	return(tmp)
}

LambdaParsing <- function(fileName){
	tallyName <- gsub(".*/","",gsub("\\.lambda.*","", fileName))
	tallyName <- gsub("_.*", "", tallyName)
	test <- tryCatch(read.table(fileName, header = F), error=function(cond){
				 message(paste(tallyName, "appears to be empty"))
				 return(data.frame(rep("k-",3), rep(NA, 3)))
})

	#test <- read.table(fileName, header = F)
	colnames(test) <- c("Metric", "Value")
	test$Sample <- tallyName
	return(test)
}

MismatchParsing <- function(fileName){
	tallyName <- gsub(".*\\/","",gsub("\\.tab","", fileName))
	tallyName <- gsub("_.*", "", tallyName)
	tmp <- read.table(fileName, header = F)
	colnames(tmp) <- c("Count", "Mismatches")
	tmp$Sample <- tallyName
	return(tmp)
		
}

DamageParsing <- function(fileName){
	tallyName <- gsub(".*\\/","",gsub("\\/(3|5).*","", fileName))
	tallyName <- gsub("_.*", "", tallyName)
	tmp <- read.delim(fileName, header = T)
	colnames(tmp) <- c("Pos", "DamageFrac")
	tmp$Sample <- tallyName	
	return(tmp)
}

FLDParsing <- function(fileName){
	tallyName <- gsub(".*\\/","",gsub("\\.tab","", fileName))
	tallyName <- gsub("_.*", "", tallyName)
	tmp <- read.delim(fileName, header = F, col.names = c("Length", "Reads"))
	tmp$Sample <- tallyName
	return(tmp)
}

NameEdits <- function(Digest){
	ifelse(Digest == "Digest1", "Digest 1",
	       ifelse(Digest == "Digest2", "Digest 2",
		      ifelse(Digest == "Digest3", "Digest 3-4",
			     ifelse(Digest == "Digest5", "Digest 5-6", "PANIC"))))
}

depurinationParsing <- function(folder){
	files <- list.files(path = folder, pattern = "lambda", full.names = T)
	Depur <- as_tibble(reduce(lapply(files, function(f){LambdaParsing(f)}), bind_rows))
	
	Depur <- Depur %>% filter(grepl("k-", Metric))
	Depur$Metric <- rep(c("ConfIntLow", "Mean", "ConfIntHigh"),4)
	Depur$Age <- 435
	Depur$Sample <- NameEdits(Depur$Sample)
	
	Depur <- Depur %>% spread(Metric, Value) #%>% left_join(age)
	#Depur[4,] <- list("Digest 1", 435, NA, NA, NA)  
	return(Depur)
}

deaminationParsing <- function(folder, age){
	files <- list.files(path = folder, pattern = "Stats_out_MCMC_iter_summ_stat.csv", full.names = T, recursive = T)
	Rate <- reduce(lapply(files, function(f){RateParsing(f)}), bind_rows)
	Rate$Age <- age
	Rate$Sample <- NameEdits(Rate$Sample)
	#Rate <- Rate %>% left_join(age)
	
	# The calculation
	Rate$Rate <- log(1/(1 - Rate$DeltaS)) * Rate$Age^-1
	
	Rate$Error <- Rate$Std/sqrt(50000) * qnorm(0.975)
	Rate$ErrorRate <- log(1/(1 - Rate$Error)) * Rate$Age^-1
	return(Rate)
}

mapDamageParsing <- function(folder){
	files <- list.files(path = folder, pattern = "5p.*", recursive = T, full.names = T)
	CT <- reduce(lapply(files, function(f){DamageParsing(f)}), bind_rows)%>% mutate(Pos = Pos - 2*Pos) %>% group_by(Sample) %>% mutate(Pos = sort(Pos))
	files <- list.files(path = folder, pattern = "3p.*", recursive = T, full.names = T)
	GA <- reduce(lapply(files, function(f){DamageParsing(f)}), bind_rows) %>% group_by(Sample) %>% mutate(Pos = sort(Pos, decreasing = T))
	
	#CT$Sample <- NameEdits(CT$Sample)
	#GA$Sample <- NameEdits(GA$Sample)

	# Trying something different
	plotDf <- CT %>% bind_rows(GA)
	return(plotDf)
}
mapDamagePlottingV2 <- function(folder, leg = F){
	files <- list.files(path = folder, pattern = "5p.*", recursive = T, full.names = T)
	CT <- reduce(lapply(files, function(f){DamageParsing(f)}), bind_rows)%>% mutate(Pos = Pos - 2*Pos) %>% group_by(Sample) %>% mutate(Pos = sort(Pos))
	files <- list.files(path = folder, pattern = "3p.*", recursive = T, full.names = T)
	GA <- reduce(lapply(files, function(f){DamageParsing(f)}), bind_rows) %>% group_by(Sample) %>% mutate(Pos = sort(Pos, decreasing = T))
	
	CT$Sample <- NameEdits(CT$Sample)
	GA$Sample <- NameEdits(GA$Sample)

#	CT <- CT %>% filter(Sample != "Digest 1")
#	GA <- GA %>% filter(Sample != "Digest 1")

	# Trying something different
	plotDf <- CT %>% bind_rows(GA)
	
	p1<- plotDf %>% ggplot(aes(x = Pos, y = DamageFrac, col = Sample)) + geom_line() + theme_bw()+
		scale_colour_manual(values = colourList) +
		xlab("Distance from Center (!?)") + ylab("Fraction Damaged") +
		coord_cartesian(ylim = c(0,0.25)) + geom_vline(xintercept = 0, lty = 2, size = 1) +
		annotate(geom = "text", x = c(-17.5,22.5), y = 0.1625, label = c("5`", "3`"), fontface = "bold")
		#scale_x_reverse()# + ggtitle("C to T Deamination")

	
	if(leg){
		p1  <- p1 + theme(legend.position = "none")
	}else{
		p1  <- p1 + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),legend.position = "none")
	}
	
	#figure <- ggdraw(figure) + draw_plot(p1sub, x = 0.125, y = 0.5, width = 0.25, height = 0.3)+ draw_plot(p2sub, x = 0.62, y = 0.5, width = 0.25, height = 0.3)
	return(p1)
}
mapDamagePlotting <- function(folder, leg = F){
	files <- list.files(path = folder, pattern = "5p.*", recursive = T, full.names = T)
	CT <- reduce(lapply(files, function(f){DamageParsing(f)}), bind_rows)
	files <- list.files(path = folder, pattern = "3p.*", recursive = T, full.names = T)
	GA <- reduce(lapply(files, function(f){DamageParsing(f)}), bind_rows)
	
	CT$Sample <- NameEdits(CT$Sample)
	GA$Sample <- NameEdits(GA$Sample)

#	CT <- CT %>% filter(Sample != "Digest 1")
#	GA <- GA %>% filter(Sample != "Digest 1")
	
	# The full plot
	p1 <- CT %>% ggplot(aes(x = Pos, y = DamageFrac, col = Sample)) + geom_line(size = 1) + theme_bw() +
		scale_colour_manual(values = colourList) + theme() +
		xlab("Distance from 5'") + ylab("Fraction Damaged") +
		coord_cartesian(ylim = c(0,0.25))# + ggtitle("C to T Deamination")
	
#	p1sub <- CT %>% ggplot(aes(x = Pos, y = DamageFrac, col = Sample)) + geom_line(size = 1) + theme_bw() +
#		scale_colour_manual(values = colour) + theme(axis.title.x = element_blank(), axis.text.y = element_blank(),
#							     axis.title.y = element_blank(), legend.position = "none") +
#		xlab("Distance from 5'") + ylab("Fraction Damaged") + coord_cartesian(xlim = c(1,5), ylim = c(0,0.35)) 
	
	
	p2 <- GA %>% ggplot(aes(x = Pos, y = DamageFrac, col = Sample)) + geom_line(size = 1) + theme_bw() +
		scale_x_reverse() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())  + coord_cartesian(ylim =c(0,0.25)) +
		scale_colour_manual(values = colourList) +
		xlab("Distance from 3'") + ylab("") # + ggtitle("G to A Deamination")
	
#	p2sub <- GA %>% ggplot(aes(x = Pos, y = DamageFrac, col = Sample)) + geom_line(size = 1) + theme_bw() +
#		scale_colour_manual(values = colour) + theme(axis.title.x = element_blank(),axis.text.y = element_blank(),
#							     axis.title.y = element_blank(), legend.position = "none") +
#		xlab("Distance from 3'") + ylab("Fraction Damaged") +
#	       	scale_x_reverse(limits = c(5,1)) + scale_y_continuous(position = "right") +
#		coord_cartesian(ylim = c(0,0.35))

	if(leg){
		p1  <- p1 + theme(legend.position = "none")
		p2  <- p2 + theme(legend.position = "none")
		figure <- ggarrange(p1,NULL,p2, nrow = 1, widths = c(1,0,1),common.legend = F, align = "hv")
	}else{
		p1  <- p1 + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),legend.position = "none")
		p2  <- p2 + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none")
		figure <- ggarrange(p1,NULL,p2, nrow = 1, widths = c(1,0,1),common.legend = F, align = "hv")
	}
	
	#figure <- ggdraw(figure) + draw_plot(p1sub, x = 0.125, y = 0.5, width = 0.25, height = 0.3)+ draw_plot(p2sub, x = 0.62, y = 0.5, width = 0.25, height = 0.3)
	return(figure)
}

sampleTrans <- function(sampleName){
	ifelse(grepl("LMN07-1", sampleName), "Digest 1",
		ifelse(grepl("LMN07-2", sampleName), "Digest 2",
		       ifelse(grepl("LMN07-3", sampleName), "Digest 3-4",
			      ifelse(grepl("LMN07-5", sampleName), "Digest 5-6",
				     ifelse(grepl("B01-3",sampleName), "Blank 3-4",
					    ifelse(grepl("B01-5", sampleName), "Blank 5-6", "Other"))))))}


#colour = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6',
#	   '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3',
#	   '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')
colour <- c("#2e294e", "#f8333c", "#007dba", "#1b998b", "#34d1bf")

colourList <- c(`Homo sapiens` = "#2e294e", `Brucella melitensis` = "#f8333c")

# Where are the files/folders that I need?
deamFolderLocation <- "NoUnmap/Ecoli/MapDamagePanGenomeV2"
depurFolderLocation <- "NoUnmap/Ecoli/Depurination"
mismatchFolderLocation <- "NoUnmap/Ecoli/MismatchesPanGenomeV2"

############# Deamination ##########
ecoli <- deaminationParsing(deamFolderLocation,435)
ecoli$Organism <- "Escherichia coli"
kaero <- deaminationParsing("NoUnmap/Kaero/MapDamageKaero", 435)
kaero$Organism <- "Klebsiella aerogenes"
human <- deaminationParsing("NoUnmap/Hsapiens/MapDamageHsapiens", 435)
human$Organism <- "Homo sapiens"

deam <- bind_rows(ecoli,human) %>% mutate(Organism =factor(Organism, levels = c("Homo sapiens","Escherichia coli","Klebsiella aerogenes"))) %>%
	ggplot(aes(x = Sample, y = Rate, ymin = Rate - ErrorRate, ymax = Rate + ErrorRate, fill = Organism)) +
	theme_bw() +scale_y_continuous(breaks = pretty_breaks(10)) +
	geom_col(position = "dodge") + geom_errorbar(colour = "black", width = 0.5, position = position_dodge(width = 0.9)) +
	scale_fill_manual(values = colour) + labs(y = "Deamination Rate") +
	theme(legend.text = element_text(face = "italic"), axis.title.x = element_blank(), axis.text.x = element_blank())
deam

ggsave(figure, file = "DeaminationPlot.png", width = 12, height = 8)

############# Depurination ##########
ecoli <- depurinationParsing(depurFolderLocation)
ecoli$Organism <- "Escherichia coli"
kaero <- depurinationParsing("NoUnmap/Kaero/Depurination")
kaero$Organism <- "Klebsiella aerogenes"
human <- depurinationParsing("NoUnmap/Hsapiens/Depurination")
human$Organism <- "Homo sapiens"

bind_rows(ecoli, kaero, human) %>% mutate(ConfIntHigh = 1/ConfIntHigh, Mean = 1/Mean, ConfIntLow = 1/ConfIntLow) %>% data.frame() 

depur <- bind_rows(ecoli, human) %>% mutate(Organism =factor(Organism, levels = c("Homo sapiens","Escherichia coli","Klebsiella aerogenes")))  %>%
	ggplot(aes(y = Mean, x = Sample, ymin = ConfIntLow, ymax = ConfIntHigh, fill = Organism)) +
	geom_col(position = "dodge") + geom_errorbar(colour = "black", width = 0.5, position = position_dodge(width = 0.9)) +
	theme_bw() + scale_fill_manual(values = colour) + scale_y_continuous(breaks = pretty_breaks(10)) +
	ylab("Depurination Rate") + theme(legend.position = "bottom")

figure <- ggarrange(ncol = 1, deam, depur, common.legend = T, legend = "bottom", align = "hv", labels = "AUTO")
ggsave(figure, file = "~/Documents/University/EcoliPaperV2/Figures/DepurandDeamPlot.png", width = 8, height = 6)

############# Readcounts ##########
files <- list.files(path = "FLD", full.names = T)
fld <- as_tibble(reduce(lapply(files, function(f){FLDParsing(f)}), bind_rows))
fld <- fld %>% mutate(Sample = replace(Sample, Sample == "BmelMapped", "Brucella melitensis"), Sample = replace(Sample, Sample == "Human", "Homo sapiens"))
fld %>% group_by(Sample,Length) %>% summarize(Reads = sum(Reads)) %>% summarize(MeanLength = mean(rep(Length, Reads)), SD = sd(rep(Length,Reads)))

MedianLength <- fld %>% group_by(Sample,Length) %>% summarize(Reads = sum(Reads)) %>% summarize(Length = median(rep(Length, Reads)))
MedianLength$Pos <- c(1e5, 3e4)
segmentPlot <- MedianLength %>% left_join(fld)


#reduce(list(ecoli,kaero, human),bind_rows) %>% mutate(Sample = NameEdits(Sample)) %>% group_by(Organism, Sample) %>%
#	summarize(MeanLength = mean(rep(Length, Reads)), StandardDeviation = sd(rep(Length, Reads)), Error = qnorm(0.975)*StandardDeviation/sqrt(sum(Reads))) %>%
#       	dplyr:::select(-c(StandardDeviation)) %>% xtable:::xtable() %>%print(file ="~/MeanLengthTable.tex")

#FLDfigureSupp <- reduce(list(kaero, human),bind_rows) %>% mutate(Sample = NameEdits(Sample)) %>%
#	mutate(Organism =factor(Organism, levels = c("Escherichia coli","Homo sapiens","Klebsiella aerogenes")))  %>%
#	ggplot(aes(x = Length, y = Reads, colour = Sample)) +
#	geom_point() + facet_grid(Organism ~.) + 
#	scale_colour_manual(values = colourList) + theme_bw() +
#	ylab("Reads") + scale_y_log10(limits = c(1,10^6), breaks = c(1,10,100,10^3,10^4, 10^5,10^6)) + scale_x_continuous(breaks = pretty_breaks(n = 10)) +
#	theme(legend.position = "bottom", strip.text.y = element_text(face = "italic")) + annotation_logticks(sides = "l") +
#	xlab("Read Length")# + ylab(bquote(log[10]("Reads")))

FLDfigure <- fld %>%
	ggplot(aes(x = Length, y = Reads, colour = Sample)) +
	geom_segment(data = segmentPlot, aes(x = Length, y = Reads, xend = Length, yend = 0), lty = 2) +
	geom_point() + #facet_grid(Organism ~.) + 
	scale_colour_manual(values = colourList) + theme_bw() +
	ylab("Reads") + scale_y_log10(limits = c(1,10^6), breaks = c(1,10,100,10^3,10^4, 10^5,10^6)) + scale_x_continuous(breaks = pretty_breaks(n = 10)) +
	theme(legend.position = "bottom", legend.text = element_text(face = "italic")) + annotation_logticks(sides = "l") +
	xlab("Read Length")  +
	annotate(geom = "text", x = 100, y = 3e5, label = "Median Read Lengths") +# + ylab(bquote(log[10]("Reads")))
	annotate(geom = "text", x = 100, y = MedianLength$Pos, label = MedianLength$Length, colour = colourList[c(2,1)])# + ylab(bquote(log[10]("Reads")))
FLDfigure


ggsave(FLDfigure, file = "FLDLog.png", width = 6, height = 4)

############# Mismatches ##########
files <- list.files("Mismatches", full.names = T)
mismatches <- as_tibble(reduce(lapply(files, function(f){MismatchParsing(f)}), bind_rows))
mismatches <- mismatches %>% mutate(Sample = replace(Sample, Sample == "BmelMapped", "Brucella melitensis"), Sample = replace(Sample, Sample == "Human", "Homo sapiens"))

figure <- mismatches %>% group_by(Sample) %>% 
        mutate(Prop = Count/sum(Count)) %>%	
	mutate(PropLabel = paste0(round(Prop * 100,3), "%")) %>%
	ggplot(aes(x = Mismatches, y = Count, color = Sample, label = PropLabel)) + geom_point() + geom_line() +
	geom_text_repel(show.legend = F, color = "black") +
	#scale_y_continuous(breaks = scales:::breaks_pretty(n = 10)) +
	scale_y_log10() + annotation_logticks(sides = "l") +
	scale_color_manual(values = colour) + theme_bw() +ylab("Mapped Reads") +
	theme(legend.text = element_text(face = "italic"), legend.position = "bottom", axis.title.y = element_blank())
figure
ggsave("BmelMismatches.png", width = 6, height = 4)
#
#figure <- reduce(list(ecoli, human),bind_rows) %>%
#	mutate(Sample = gsub("Mimatches|Mismatches","", Sample), Organism = factor(Organism, levels = c("Homo sapiens", "Escherichia coli"))) %>%
#	mutate(Sample = NameEdits(Sample)) %>% 
#	ggplot(aes(x = Mismatches, y = Count, fill = Sample)) + geom_col() +
#	facet_grid(Organism ~ Sample) +
#	scale_y_log10() + annotation_logticks(sides = "l") +
#	scale_fill_manual(values = colour) + theme_bw() +
#	theme(strip.text.y = element_text(face = "italic")) +
#	ylab("Reads") #+ scale_y_continuous(breaks = pretty_breaks(n = 10)) 

############# Overlapping Smiles ##########
#ecoli <- mapDamagePlottingV2(deamFolderLocation,leg = F)
mapDamage <- mapDamageParsing("MapDamage")
mapDamage <- mapDamage %>% mutate(Sample = replace(Sample, Sample == "BrucellaMelitensis", "Brucella melitensis"), Sample = replace(Sample, Sample == "HomoSapiens", "Homo sapiens"))

annotationDF <- mapDamage %>% group_by(Sample) %>% filter(Pos == -25)
annotationDF$Pos <- c(0.15,0.175)
#annotationDF$Pos <- c(0.1,0.125,0.15,0.175)
annotationDF <- annotationDF %>% mutate(Sample =factor(Sample, levels = c("Homo sapiens","Brucella melitensis")))


mapDamagefigureSupp <- reduce(list(kaero, human),bind_rows) %>% filter(Sample != "Digest 1") %>%
	mutate(Organism =factor(Organism, levels = c("Escherichia coli","Homo sapiens","Klebsiella aerogenes"))) %>%
	ggplot(aes(x = Pos, y = DamageFrac, col = Sample)) + geom_line() + theme_bw()+
	scale_colour_manual(values = colourList) + facet_grid(Organism ~.) +
	xlab("Distance from Center") + ylab("Fraction Damaged") +
	coord_cartesian(ylim = c(0,0.25)) + geom_vline(xintercept = 0, lty = 2, size = 1) +
	theme(strip.text.y = element_text(face = "italic"))

mapDamagefigure <- mapDamage  %>%
	mutate(Sample =factor(Sample, levels = c("Homo sapiens","Brucella melitensis"))) %>%
	ggplot(aes(x = Pos, y = DamageFrac, col = Sample)) + geom_line() + theme_bw()+
	scale_colour_manual(values = colourList) +# facet_grid(Organism ~.) +
	xlab("Distance from Center") + ylab("Fraction Damaged") +
	coord_cartesian(ylim = c(0,0.25)) + geom_vline(xintercept = 0, lty = 2, size = 1) +
	theme(legend.text = element_text(face = "italic"), legend.position = "none") +
	annotate(geom = "text", x = -10, y = 0.20, label = "Maximum CT Deamination") +
	annotate(geom = "text", x = -10, y = annotationDF$Pos, label = signif(annotationDF$DamageFrac,3), colour = colourList[c(2,1)])

ggsave("MapDamageFigure.pdf", width = 6, height = 4)
	#annotate(geom = "text", x = 10, y = , label = c("5`", "3`"), fontface = "bold")

bottom <- ggarrange(FLDfigure, figure, common.legend =T, legend = "bottom")
ggarrange(mapDamagefigure, bottom, ncol = 1)
ggsave("DamagePlots.pdf", width = 12, height = 8)
