library(STRINGdb)
library(dplyr)

string_db <- STRINGdb$new(version="11.5", score_threshold = 700, species = 224914, input_directory = "")

# Now to get the enrichment results
genes <- read.table("CoreMissingAllGeridu.list", header =F, col.names = "Gene")
mappedGenes <- string_db$map(genes, "Gene", removeUnmappedRows = T)
GeriduenrichedFunctions <- string_db$get_enrichment(mappedGenes$STRING_id) %>% as_tibble()

# Now to compare the two sets and see what's unique!
corsiniEnriched <- CorsinEenrichedFunctions %>% select(`category`, term, description)
geriduEnriched <- GeriduenrichedFunctions %>% select(`category`, term, description)

differencesCor <- setdiff(corsiniEnriched, geriduEnriched)
differencesCor$Sample <- "Brancorsini"
differencesGer <- setdiff(geriduEnriched, corsiniEnriched)
differencesGer$Sample <- "Geridu"

#
differencesCor %>% bind_rows(differencesGer) %>% write.table(file = "DifferentiallyExpressedFunctions.tab", sep = "\t", quote =F, row.names = F)
