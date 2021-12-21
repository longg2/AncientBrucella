library(ggplot2)
library(ggpubr)
library(ggrepel)
library(dplyr)
library(tidyr)
library(cluster)

tmp <- read.table("SnpDistances.tab", header =T) 
colnames(tmp) <- gsub("^X", "", colnames(tmp))
colnames(tmp) <- gsub("\\.", "-", colnames(tmp))

snpDist <- dist(tmp)

fit <- cmdscale(snpDist, eig =T, k = 4, add = T)
coord <- fit$points %>% as_tibble()
coord$Genome <- rownames(fit$points)
contrib <- fit$eig/sum(fit$eig) * 100
coord$Colour <- ifelse(coord$Genome == "MappedReads", "Ancient", "Modern")

p1 <- coord %>%
       	ggplot(aes(x = V1, y = V2, label = Genome, colour = Colour)) +
	geom_hline(yintercept=0, lty = 2, colour = "grey90") + 
	geom_vline(xintercept=0, lty = 2, colour = "grey90") +
	geom_point() +
	xlab(bquote("PCoA 1 ("~.(round(contrib[1],2))~"%)")) +
	ylab(bquote("PCoA 2 ("~.(round(contrib[2],2))~"%)")) +
	scale_colour_manual(values = c("red", "black")) +
	geom_text_repel(show.legend = F) +
	theme_classic()
p1

p2 <- coord %>%
       	ggplot(aes(x = V3, y = V4, label = Genome, colour = Colour)) +
	geom_hline(yintercept=0, lty = 2, colour = "grey90") + 
	geom_vline(xintercept=0, lty = 2, colour = "grey90") +
	geom_point() +
	xlab(bquote("PCoA 3 ("~.(round(contrib[3],2))~"%)")) +
	ylab(bquote("PCoA 4 ("~.(round(contrib[4],2))~"%)")) +
	scale_colour_manual(values = c("red", "black")) +
	geom_text_repel(show.legend = F) +
	theme_classic()
p2

ggarrange(p1,p2, legend = "bottom", align = "hv", ncol = 1, common.legend = T)

ggsave("SNPPCoA.pdf")
