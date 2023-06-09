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
library(latex2exp)

pvalText <- function(pval){

	if(pval < 0.001){
		return(paste0("$p < ", 0.001, "$"))
	}else{
		return(paste0("$p = ", round(pval, 3), "$"))
	}
}

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
mapDamageParsingOld <- function(folder){
	files <- list.files(path = folder, pattern = "5p.*", recursive = T, full.names = T)
	#CT <- reduce(lapply(files, function(f){DamageParsing(f)}), bind_rows)%>% mutate(Pos = Pos - 2*Pos) %>% group_by(Sample) %>% mutate(Pos = sort(Pos))
	CT <- reduce(lapply(files, function(f){DamageParsing(f)}), bind_rows)%>% group_by(Sample) %>% mutate(Deam = "5'")
	files <- list.files(path = folder, pattern = "3p.*", recursive = T, full.names = T)
	GA <- reduce(lapply(files, function(f){DamageParsing(f)}), bind_rows) %>% group_by(Sample) %>% mutate(Deam = "3'")
	
	# Trying something different
	plotDf <- CT %>% bind_rows(GA)
	return(plotDf)
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

colour = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6',
	   '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3',
	   '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')
#colour <- c("#2e294e", "#f8333c", "#007dba", "#1b998b", "#34d1bf")

colourList <- c("Human" = colour[4], "Human Trimmed" = colour[7], "Brucella melitensis" = colour[1])
############# Readcounts ##########
files <- list.files(path = "FLD", full.names = T)
fld <- as_tibble(reduce(lapply(files, function(f){FLDParsing(f)}), bind_rows))
fld <- fld %>% filter(!grepl("Nodule1", Sample)) %>% 
	mutate(Sample = ifelse(Sample == "Brancorsini", "Brucella melitensis", ifelse(Sample == "HumanBrancorsini","Human", "Human Trimmed")))
fld %>% group_by(Sample,Length) %>% summarize(Reads = sum(Reads)) %>% summarize(MeanLength = mean(rep(Length, Reads)), SD = sd(rep(Length,Reads)))

MedianLength <- fld %>% group_by(Sample,Length) %>% summarize(Reads = sum(Reads)) %>% summarize(Length = median(rep(Length, Reads)))
MedianLength$Pos <- c(1e5, 3e4, 1e4)
segmentPlot <- MedianLength %>% left_join(fld)

# Quick test to see if I can get a significant difference
fldModel <- fld %>% filter(!grepl("Geridu", Sample)) %>% lm(formula = log10(Reads) ~ Length * Sample)
par(mfrow = c(2,2))
plot(fldModel)
dev.off()
fldModel %>% summary()

rsquaredAdjust <- summary(fldModel)$adj.r.squared
f <- summary(fldModel)$fstatistic
pval <- pf(f[1], f[2], f[3], lower.tail = F)
pvalLab <- pvalText(pval)

BranSlope <- abs(coef(fldModel)["Length"]) # We can rearrange the formula so that it works out
HumanSlope <- abs(sum(coef(fldModel)[c("Length:SampleHuman","Length")]))
HumanTrimmedSlope <- abs(sum(coef(fldModel)[c("Length:SampleHuman Trimmed","Length")]))

FLDfigure <- fld %>%
	ggplot(aes(x = Length, y = Reads, colour = Sample)) +
	geom_segment(data = segmentPlot, aes(x = Length, y = Reads, xend = Length, yend = 0), lty = 2, show.legend = F) +
	geom_point() + #facet_grid(Organism ~.) + 
#	geom_smooth(method = "rlm") +
	geom_text(data = MedianLength, aes(x = 90, y = Pos, label = Length, colour = Sample),show.legend = F) +# + ylab(bquote(log[10]("Reads")))
	scale_colour_manual(values = colourList) + theme_classic() +
	ylab("Reads") + scale_y_log10(limits = c(1,10^6), breaks = c(1,10,100,10^3,10^4, 10^5,10^6)) + scale_x_continuous(breaks = pretty_breaks(n = 10)) +
	theme(legend.text = element_text(face = "italic"), legend.position = "bottom") +
	annotation_logticks(sides = "l") +
	xlab("Read Length")  +
	annotate(geom = "text", x = 90, y = 3e5, label = "Median Read Lengths") +
	annotate(geom = "text", x = 60, y = 1e2, label = TeX(paste0("$R^2_{adj} = ", round(rsquaredAdjust, 3), "$"))) +
	annotate(geom = "text", x = 60, y = 3e1, label = TeX(pvalLab)) +
	annotate(geom = "text", x = 120, y = 1e4, label = "Depurination") +
	annotate(geom = "text", x = 120, y = 3e3, label = round(BranSlope, 3), colour = colour[1]) +
	annotate(geom = "text", x = 120, y = 1e3, label = round(HumanSlope, 3), colour = colour[4]) +
	annotate(geom = "text", x = 120, y = 3e2, label = round(HumanTrimmedSlope, 3), colour = colour[7]) 
	#guides(colour = guide_legend(nrow = 2))# + ylab(bquote(log[10]("Reads")))
FLDfigure
ggsave(FLDfigure, file = "FLDLogSupplement.pdf", width = 6, height = 4)
#ggsave(FLDfigure, file = "FLDLog.pdf", width = 6, height = 4)

############# Mismatches ##########
files <- list.files("Mismatches", full.names = T)
mismatches <- as_tibble(reduce(lapply(files, function(f){MismatchParsing(f)}), bind_rows))
mismatches <- mismatches %>% filter(!grepl("Nodule1", Sample)) %>% 
	mutate(Sample = ifelse(Sample == "Brancorsini", "Brucella melitensis", ifelse(Sample == "HumanBrancorsini","Human", "Human Trimmed")))
#mismatches <- mismatches %>% 
#	mutate(Sample = replace(Sample, Sample == "Brancorsini", "Brucella melitensis - Brancorsini"), Sample = replace(Sample, Sample == "Nodule1", "Brucella melitensis - Geridu"), Sample = replace(Sample, Sample == "HumanBrancorsini", "Human - Brancorsini"), Sample = replace(Sample, Sample == "HumanNodule1", "Human - Geridu")) %>%
#	mutate(TMP = ifelse(grepl("Brancorsini", Sample), "Brancorsini", "Geridu"))
#mismatches <- mismatches %>% mutate(Sample = replace(Sample, Sample == "BmelMapped", "Brucella melitensis"), Sample = replace(Sample, Sample == "Human", "Homo sapiens"))

mismatchMeans <- mismatches %>% group_by(Sample) %>% summarize(MeanMismatches = sum(Count * Mismatches)/sum(Count), MedianMismatches = median(rep(Mismatches, Count))) %>% mutate(height = c(1e6,3e5, 1e5))

figure <- mismatches %>% group_by(Sample) %>% 
	ggplot(aes(x = Mismatches, y = Count, fill = Sample)) + geom_col(position = "dodge") +
	geom_text(data = mismatchMeans, aes(x = 4, y = height, label = round(MeanMismatches,3), colour = Sample),show.legend = F) +# + ylab(bquote(log[10]("Reads")))
	scale_x_continuous(breaks = scales:::breaks_pretty(n = 9)) +
	scale_y_log10() + annotation_logticks(sides = "l") +
	scale_colour_manual(values = colourList) + scale_fill_manual(values = colourList) +
	theme_classic() + ylab("Mapped Reads") +
	theme(legend.text = element_text(face = "italic"), legend.position = "bottom") +
	annotate(geom = "text", x = 4, y = 3e6, label = "Mean Mismatches") +
	ylab("Reads") 
	#guides(colour = guide_legend(nrow = 2))# + ylab(bquote(log[10]("Reads")))
figure
ggsave("BmelMismatches.png", width = 9, height = 6)
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
mapDamage <- mapDamageParsingOld("MapDamage") 
mapDamage <- mapDamage %>% filter(!grepl("Nodule1", Sample)) %>% 
	mutate(Sample = ifelse(Sample == "Brancorsini", "Brucella melitensis", ifelse(Sample == "HumanBrancorsini","Human", "Human Trimmed")))
#mapDamage <- mapDamage %>% mutate(Sample = replace(Sample, Sample == "Brancorsini", "Brucella melitensis - Brancorsini"), Sample = replace(Sample, Sample == "Nodule1", "Brucella melitensis - Geridu"), Sample = replace(Sample, Sample == "HumanBrancorsini", "Human - Brancorsini"), Sample = replace(Sample, Sample == "HumanNodule1", "Human - Geridu"))

annotationDF <- mapDamage %>% group_by(Sample) %>% filter(Pos == 1)
annotationDF$Pos <- rep(c(0.15,0.175),2)

fivePrime <- mapDamage %>% filter(Deam == "5'") %>%
	#mutate(Sample =factor(Sample, levels = c("Homo sapiens","Brucella melitensis"))) %>%
	ggplot(aes(x = Pos, y = DamageFrac, col = Sample)) + geom_line() + theme_classic()+
	scale_colour_manual(values = colourList) +# facet_grid(Organism ~.) +
	xlab("Distance from End") + ylab("Fraction Damaged") +
	#geom_text(data = annotationDF, aes(x = c(rep(10,8)), y = Pos, label = signif(DamageFrac, 3), colour = Sample),show.legend = F) +# + ylab(bquote(log[10]("Reads")))
	geom_text(data = annotationDF %>% filter(Deam == "5'"), aes(y = Pos, x = 12.5, label = signif(DamageFrac, 3), colour = Sample), show.legend = F) +
	coord_cartesian(ylim = c(0,0.25)) +
	theme(legend.position = "bottom") +
	annotate(geom = "text", x = 12.5, y = 0.2, label = "CT Deamination")

threePrime <- mapDamage %>% filter(Deam != "5'") %>%
	#mutate(Sample =factor(Sample, levels = c("Homo sapiens","Brucella melitensis"))) %>%
	ggplot(aes(x = Pos, y = DamageFrac, col = Sample)) + geom_line() + theme_classic()+
	scale_colour_manual(values = colourList) +# facet_grid(Organism ~.) +
	xlab("Distance from End") + ylab("Fraction Damaged") +
	scale_x_reverse() + scale_y_continuous(position = "right")+
	geom_text(data = annotationDF %>% filter(Deam != "5'"), aes(y = Pos, x = 12.5, label = signif(DamageFrac, 3), colour = Sample), show.legend = F) +
	coord_cartesian(ylim = c(0,0.25)) +
	theme(legend.position = "bottom", axis.title.y = element_blank(), axis.text.y = element_blank()) +
	annotate(geom = "text", x = 12.5, y = 0.2, label = "GA Deamination")
#mapDamagefigure <- ggarrange(fivePrime, threePrime, common.legend = T, legend = "bottom", labels = c("5`","3`"), align = "hv")
mapDamagefigure <- ggarrange(fivePrime, threePrime, common.legend = T, legend = "none", align = "hv")
ggsave("MapDamageFigure.pdf", width = 6, height = 4)
	#annotate(geom = "text", x = 10, y = , label = c("5`", "3`"), fontface = "bold")

# Committee Meeting
bottom <- ggarrange(FLDfigure, figure, common.legend =T,align = "hv", legend = "bottom", labels = c("b", "c"))
ggarrange(mapDamagefigure, bottom, labels = c("a",""), nrow = 2)
ggsave("Figure1CM.pdf", width = 9, height = 6)
