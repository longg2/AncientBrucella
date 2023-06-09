#! /usr/bin/Rscript 
library(dplyr)
library(tidyr)
library(ggplot2)
library(pbapply)
#kneedle <- function(x,guess=length(x)){
#       	# This is an implementation by Zachery Dickson. Comes from https://raghavan.usc.edu/papers/kneedle-simplex11.pdf 
#    if(is.na(guess) || guess > length(x) / 2){
#        guess = length(x)/2
#    }
#    x = x[1:(guess*2)]
#    y = lowess(x)$y
#    y0 = y[1]
#    yn= y[length(y)]
#    x0=0
#    xn=length(y)-1
#    m = (yn-y0)/(xn-x0)
#    b = y0
#    d = abs(y - (m*(x0:xn) + b)) * asin(pi/4)
#    return(which.max(d) + 1)
#}

#fileName <- "SubspeciesStrict.log"
arguments <- commandArgs(trailingOnly = T)
fileName <- arguments[1]

if(!grepl("\\.log", fileName)){
	stop("This isn't a BEAST log file")
}else if(!file.exists(gsub("\\.log.*","\\.trees", fileName))){
	stop("Missing the trees!")
}
ESSthresh <- 200

logFile <- read.delim(comment.char = "#", file = fileName)
#logList <- logFile[,-1] %>% as.list()
#%>% as_tibble()

stateSize <- logFile$state[3] - logFile$state[2]
guessSize <- logFile$state[nrow(logFile)] * 0.2
#guessSize <- logFile$state[nrow(logFile)] * 0.2 / stateSize

#burnInPeriods <- pblapply(logList, cl = 4, function(x){kneedle(x, guess = guessSize) * stateSize})
#burnInPeriods <- pbapply(logFile %>% filter(state > 1) %>% select(-state), MARGIN = 2, function(x){kneedle(x, guess = guessSize) * stateSize})

tmp <- system(paste("loganalyser -burnin", formatC(guessSize, format = "d"),fileName), ignore.stderr = T, intern = T)[-c(1:17)]
#tmp <- system(paste("loganalyser -burnin",burnInPeriods["joint"],fileName), ignore.stderr = T, intern = T)[-c(1:17)]
tmp <- gsub("\t\\*|\t$","",tmp)
beg <- which(nchar(tmp) == 0)
logAnalyserOut <- tmp[-c(beg:length(tmp))] %>% strsplit(split = "\t") %>% list2DF() %>% t() %>% as_tibble()

colnames(logAnalyserOut) <- logAnalyserOut[1,]
logAnalyserOut <- logAnalyserOut[-1,] %>% mutate(across(!statistic, as.numeric))

# Let's plot all of these results
#burnInPlot <- tibble("Statistic" = names(burnInPeriods), "state" = unname(burnInPeriods))
#plotData <- logFile %>% select(state, joint, prior, likelihood, pInv, alpha, age.root., treeLikelihood) %>% 
plotData <- logFile %>% select(state, joint, prior, likelihood, pInv, alpha, ucld.mean,ucld.stdev,meanRate, age.root., treeLikelihood) %>% 
#plotData <- logFile %>% select(state, joint, prior, likelihood, pInv, alpha, meanRate, age.root., treeLikelihood) %>% 
	filter(state > 10000) %>% pivot_longer(-state, names_to = "Statistic", values_to = "Value")
traceplotData <- plotData %>% filter(state %% 100000 == 0)

tracePlot <- traceplotData %>% ggplot(aes(x = state, y = Value)) +
	geom_line() + theme_classic() + facet_wrap("Statistic", scales = "free_y") + ggtitle("Traces filtering out the first 10000 states") +
	geom_smooth() +
	#annotate(geom = "vline", xintercept = burnInPeriods["joint"], colour = "red", lty =2)
	geom_vline(xintercept = guessSize, colour = "red", lty = 2) +
	theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggsave(tracePlot, file = gsub(".log","TracePlot.png", fileName), width = 24, height = 18)

summaryLog <- logFile %>%
       	#select(state, joint, prior, likelihood, pInv, alpha, meanRate, age.root., treeLikelihood) %>% 
#       	select(state, joint, prior, likelihood, pInv, alpha, age.root., treeLikelihood) %>% 
       	select(state, joint, prior, likelihood, pInv, alpha, ucld.mean,ucld.stdev,meanRate, age.root., treeLikelihood) %>% 
	pivot_longer(-state, names_to = "Statistic", values_to = "Value") %>%
       	mutate(state = replace(state, state < guessSize,NA)) %>%
       	#mutate(state = replace(state, state < burnInPeriods["joint"],NA)) %>%
       	filter(!is.na(state)) %>% group_by(Statistic) %>%
       	summarize(Mean = mean(Value), Median = median(Value), SD = sd(Value), SE = qnorm(0.975) * sd(Value)/sqrt(length(Value)),
	CIHi = Mean + SE, CILo = Mean - SE)

probPlot <- plotData %>%
       	mutate(state = replace(state, state < guessSize,NA)) %>%
       	#mutate(state = replace(state, state < burnInPeriods[["joint"]],NA)) %>%
       	filter(!is.na(state)) %>% ggplot(aes(x = Value)) +
	geom_density(colour = "black", fill = "lightblue") + 
	geom_rect(data = summaryLog, aes(ymin = -Inf, ymax = Inf, xmin = Mean - 1.96 * SD, xmax = Mean + 1.96 * SD), inherit.aes = F, alpha = 0.3) +
	geom_vline(data = summaryLog, aes(xintercept = Mean), colour = "red", lty = 2) +
	geom_vline(data = summaryLog, aes(xintercept = Median), colour = "blue", lty = 2) +
	#geom_rect(data = summaryLog, inherit.aes = F, aes(xmin = CILo, xmax = CIHi, ymin = -Inf, ymax = Inf), colour = "grey90", alpha = 0.1) +
	theme_classic() + facet_wrap("Statistic", scales = "free") + 
	ggtitle(paste0("Posterior Densities after Burn-In = ", formatC(guessSize, format = "d"))) +
	#ggtitle(paste0("Posterior Densities after Burn-In = ", burnInPeriods["joint"])) +
	theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
	#geom_vline(data = burnInPlot, aes(xintercept = state), colour = "red", lty = 2)
ggsave(probPlot, file = gsub(".log","ProbPlot.png", fileName), width = 24, height = 18)

# Need to report ESSs which are lower than 200. Depending on how bad it is, we want to stop right here so that we don't make a useless tree
notSampledEnough <- logAnalyserOut %>% filter(ESS < ESSthresh)
if(nrow(notSampledEnough) > 0){
	write.table(notSampledEnough, file = paste0(gsub(".log","",fileName),"ESSLow.log"), sep = "\t", quote = F, row.names = F)
	if("joint" %in% notSampledEnough$statistic){
		stop("Joint not sampled enough!")
	}
}

# Let's get that tree made
system(paste("treeannotator -burnin",formatC(guessSize, format = "d"), "-heights mean",gsub("log","trees",fileName),gsub("log","nexus",fileName)))
#system(paste("treeannotator -burnin",burnInPeriods["joint"], "-heights mean",gsub("log","trees",fileName),gsub("log","nexus",fileName)))
#system(paste("treeannotator -burnin",burnInPeriods["joint"], "-heights ca",gsub("log","trees",fileName),gsub("log","nexus",fileName)))

