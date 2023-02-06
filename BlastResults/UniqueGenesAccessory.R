library(dplyr)
library(tidyr)

# Reading in the Accesory Gene Functions
AccessFunctions <- read.delim("AccessoryGenesUnique.tab") %>% as_tibble() %>% select(-c(PMID,PDBID,COGMembershipClass)) %>% distinct()

# Reading in the Gene Cluster P/A
clusterPA <- read.delim("../AssemblyOnly/AccessoryGeneLocation.tab") %>% as_tibble()

# Merging the results together
Alltogether <- clusterPA %>% full_join(AccessFunctions, by = c("Gene" = "Query")) %>% filter(!is.na(Africa.America))

# Now, we're only interested in those accessory genes that aren't present everywhere
differential <- Alltogether %>% filter(!(Africa.America & Fertile.Crescent & Indo.Pacific & Russia & Western.Mediterranean))

write.table(differential, row.names = F, quote = F, file = "AccessoryGenesUniqueFunctions.tab", sep = "\t")
