library(dplyr)
library(parallel)
library(pbapply)

# Reading the results
op <- pboptions(type = "timer")
ncores <- 7

classifiedHits <- read.delim(paste0("IdentifiedProteinsAfrica.list"), header = T)
	
for (file in list.files("WeirdAfrica/", full.names = T)){
	#FileName
	filename <- gsub(".tab*", "",basename(file))
	tmp <- read.delim(file, header = F,
			  col.names = c("Query", "Match", "PIdent", "Length", "Mismatches", "Gapopen", "QStart", "QEnd", "SStart", "SSend", "Evalue", "Bitscore", "Taxa")) %>% as_tibble()
	
	tmpList <- split(tmp, as.factor((tmp$Query)))
	
	FilteredAmbiguousGenes<- pblapply(cl = 8,tmpList, function(x){ 
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
			return(filtResults)
	 }) %>% bind_rows() %>% filter(PIdent >= 90)

	FilteredAmbiguousGenes <- FilteredAmbiguousGenes %>% mutate(Match = gsub("\\.[0-9]$","", Match)) %>% left_join(classifiedHits, c("Match" = "ProteinID")) %>%
		select(Query, Match, Name, Organism) %>% distinct()
	#write.table(FilteredAmbiguousGenes,file = paste0("COGCoreFixed/",filename,".tab"),row.name = F, quote = F, sep = "\t")
	
	#FilteredAmbiguousGenes %>% pull(Match) %>% unique() %>% write.table(paste0("BlastTop/",filename,".list"), col.names = F, quote = F, row.names = F)
}
	write.table(FilteredAmbiguousGenes, file = "WeirdAfrica.tab", sep = "\t", quote = F, row.name = F)
