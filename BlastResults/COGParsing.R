library(dplyr)
library(parallel)
library(pbapply)

# Reading the COG data
prot2cog <- read.csv("cog-20.cog.csv.gz", header =F,
	 col.names= c("Gene", "Assembly", "Protein", "ProteinLength",
		      "COGFootprintCoord", "COGFootprintLength", "COG",
		      "reserved","COGMembershipClass","Bitscore", "E-Value",
		      "COGProfileLength", "ProteinFootprintCoord")) %>% as_tibble() %>%
	mutate(Protein = gsub("\\.[0-9]$","",Protein))

cog2fun <-  read.delim("cog-20.def.tab", header =F,
	 col.names= c("COG", "FunctionalCategory", "Name", "Gene",
		      "FunctionalPathway", "PMID", "PDBID")) %>% as_tibble()


# Reading the results
op <- pboptions(type = "timer")
ncores <- 6
	
for(file in list.files("COGInputs", full.names = T)){
	#FileName
	filename <- gsub(".tab.gz", "",basename(file))
	tmp <- read.delim(file, header = F,
			  col.names = c("Query", "Match", "PIdent", "Length", "Mismatches", "Gapopen", "QStart", "QEnd", "SStart", "SSend", "Evalue", "Bitscore", "Taxa")) %>% as_tibble()
	
	tmpList <- split(tmp, as.factor((tmp$Query)))
	
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
	 }) %>% bind_rows() %>% filter(PIdent >= 90) %>% mutate(Match = gsub("_[0-9]$","",Match))

	FilteredAmbiguousGenes <- FilteredAmbiguousGenes %>% select(Query,Match) %>% left_join(prot2cog %>% select(Protein,COG,COGMembershipClass), by = c("Match" = "Protein")) %>% left_join(cog2fun) %>% distinct()

	write.table(FilteredAmbiguousGenes,file = paste0(filename,".tab"),row.name = F, quote = F, sep = "\t")
}
