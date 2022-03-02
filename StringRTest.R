# Converting to analyzing COG data
library(dplyr)
library(ggplot2)
library(ggpubr)

# The functions translated
cogFun <- read.delim("BlastResults/fun-20.tab", header = F, col.names = c("FunctionID", "Colour", "Name")) %>% as_tibble() %>% mutate(Colour = paste0("#",Colour))

# COG Colours
colour <- as.list(cogFun$Colour)
names(colour) <- cogFun$Name

# Shared core genes missing
coreAll <- read.delim("BlastResults/COGClassified/CoreGenesMissingAll.tab", header = T, na.string = "") %>% as_tibble() %>%
       	select(-c(COGMembershipClass)) %>% distinct() %>% arrange(COG) %>% mutate(PMID = as.character(PMID))

# Core genes missing in only Bran
coreCor <- read.delim("BlastResults/COGClassified/CoreGenesMissingCorsini.tab", header = T, na.string = "") %>% as_tibble() %>%
       	select(-c(COGMembershipClass)) %>% distinct() %>% arrange(COG) %>% mutate(PMID = as.character(PMID))

# Core genes missing in only Geridu
coreGer <- read.delim("BlastResults/COGClassified/CoreGenesMissingGeridu.tab", header = T, na.string = "") %>% as_tibble() %>%
       	select(-c(COGMembershipClass)) %>% distinct() %>% arrange(COG) %>% mutate(PMID = as.character(PMID))

# The preliminary analysis (can we pull anything more interesting this way?)
coreCor <- coreAll %>% bind_rows(coreCor)
CorsiniFunctionCategories <- coreCor %>% pull(FunctionalCategory) %>% strsplit(split = "") %>% unlist() %>% table() %>% as_tibble()
colnames(CorsiniFunctionCategories)  <- c("FunctionID", "Count")

# Want to sort these by number of genes
CorsiniFunctionCategories <- CorsiniFunctionCategories %>% arrange(-Count) %>% left_join(cogFun %>% select(-Colour)) %>%
       	mutate(Name = factor(x = Name, levels = Name))

corsini <- ggplot(CorsiniFunctionCategories, aes(x = Name, fill = Name, y = Count)) + geom_col() +
	scale_fill_manual(values = colour) + theme_classic() +
	theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank()) +
	guides(fill = guide_legend(ncol = 1, title = "COG Category"))

# The preliminary analysis (can we pull anything more interesting this way?)
coreGer <- coreAll %>% bind_rows(coreGer)
GeriduFunctionCategories <- coreGer %>% pull(FunctionalCategory) %>% strsplit(split = "") %>% unlist() %>% table() %>% as_tibble()
colnames(GeriduFunctionCategories)  <- c("FunctionID", "Count")

# Want to sort these by number of genes
GeriduFunctionCategories <- GeriduFunctionCategories %>% arrange(-Count) %>% left_join(cogFun %>% select(-Colour)) %>%
       	mutate(Name = factor(x = Name, levels = Name))

geridu <- ggplot(GeriduFunctionCategories, aes(x = Name, fill = Name, y = Count)) + geom_col() +
	scale_fill_manual(values = colour) + theme_classic() +
	theme(axis.text.x = element_blank()) +
	guides(fill = guide_legend(ncol = 1))

ggarrange(corsini, geridu, labels = "AUTO", common.legend = T, ncol = 1, legend = "right")

ggsave("COGFunctions.png", width = 12, height = 9)
