library(dplyr)
library(parallel)
library(pbapply)

# Reading the results
op <- pboptions(type = "timer")
ncores <- 10
tmp <- list.files("BlastOut", full.names = T)
	
for (file in list.files("BlastOut", full.names = T)){
	#FileName
	filename <- gsub(".tab", "",basename(file))
	tmp <- read.delim(file, header = F,
			  col.names = c("Query", "Match", "PIdent", "Length", "Mismatches", "Gapopen", "QStart", "QEnd", "SStart", "SSend", "Evalue", "Bitscore", "Taxa")) %>% as_tibble()
	
	tmpList <- split(tmp, as.factor((tmp$Query)))
	
	FilteredAmbiguousGenes<- pblapply(cl = 10,tmpList, function(x){ 
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
	
	classifiedHits <- read.delim(paste0("ProteinFiguredOut/",filename,".tab"), header = F, col.names = c("Match", "Name", "Organism"))
	FilteredAmbiguousGenes <- FilteredAmbiguousGenes %>% mutate(Match = gsub("\\.[0-9]$","", Match)) %>% left_join(classifiedHits) %>%
		select(Query, Match, Name, Organism)
	write.table(FilteredAmbiguousGenes,file = paste0("Classified/",filename,".tab"),row.name = F, quote = F, sep = "\t")

	
	#FilteredAmbiguousGenes %>% pull(Match) %>% unique() %>% write.table(paste0("BlastTop/",filename,".list"), col.names = F, quote = F, row.names = F)
}


# Assuming its been done
#classifiedHits <- read.delim("IdentifiedGenes.list", header = F, col.names = c("Match", "Name", "Organism"))
#FilteredAmbiguousGenes <- FilteredAmbiguousGenes %>% mutate(Match = gsub("\\.[0-9]$","", Match))
#tmp <- FilteredAmbiguousGenes %>% left_join(classifiedHits) %>% filter(grepl(" VI |T6SS|type VI", Name))
#FilteredAmbiguousGenes %>% left_join(classifiedHits) %>% select(Query, Match, Name, Organism) %>% 
#	write.table(file = "ClassifiedGenes.tab",row.name = F, quote = F, sep = "\t")
#
#tmp %>% pull(Query) %>% write.table(file = "T6SSGenesFound.list", col.name = F, row.name = F, quote = F)

