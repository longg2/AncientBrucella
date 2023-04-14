library(NELSI)
library(lubridate)
library(phytools)

trName <- 'AssemblyMLTrees/GlobalReview.boottrees'
datesName <- '../BruceDatesUpdated.tab'

tr1 <- read.tree(trName)
dates1 <- read.table(datesName, head = F, stringsAsFactors = F)

dates1$newname <- NA
head(dates1)
for(i in 1:nrow(dates1)){
    if(nchar(dates1$V2[i]) ==  4){
        dates1$V2[i] <- paste0(dates1$V2[i], '-06-15')
    }else if(nchar(dates1$V2[i]) ==  7 | nchar(dates1$V2[i]) == 6){
        dates1$V2[i] <- paste0(dates1$V2[i], '-15')
    }
    dates1$newname[i] <- paste(dates1$V1[i],
                               round(decimal_date(ymd(dates1$V2[i])), 3),
                               sep = '_')
}
dates1

#outgroupLabel <- dates1$newname[grepl("Reference",dates1$V1)] # This is for WMed
outgroupLabel <- dates1$newname[grepl("abortus",dates1$V1)] # This is for the Global

datedTrees <- list()
rates <- vector()
tMRCA <- vector()
obFunc <- vector()
for(i in 1:length(tr1)){
    trTemp <- tr1[[i]]
    trTemp$tip.label <- gsub('[.].+', '', trTemp$tip.label)
    trTemp$tip.label <- dates1$newname[match(trTemp$tip.label, dates1$V1)]
    trTemp <- root(trTemp, outgroup = outgroupLabel, resolve.root = T)
    write.tree(trTemp, 'temp.wdates.tree')
    cat('replicate', i, '\n')
    make.lsd.dates(trTemp, outfile = 'outdateTemp.date')
    system(paste('~/Applications/lsd2/bin/lsd2_unix -i temp.wdates.tree -d outdateTemp.date -r l -s 1500000 -o resultTemp'))
    rates[i] <- as.numeric(gsub(',.+', '',
                                gsub('.+rate ', '',
                                     grep(' rate.+tMRCA', readLines('resultTemp'), value = T))))
    tMRCA[i] <- as.numeric(gsub(' ,.*', '',
                                gsub('.+tMRCA ', '',
                                     grep('tMRCA.+,', readLines('resultTemp'), value = T))))
    obFunc[i] <- as.numeric(gsub('.*function ', '',
                                     grep('tMRCA.+,', readLines('resultTemp'), value = T)))
    if(is.na(tMRCA[i])) stop("Check")
    datedTrees[[i]] <- read.nexus('resultTemp.date.nexus')
}
write.nexus(datedTrees, file = 'datedBootAssemblyGlobalReview.nexus')
system("rm -rf tmp.wdates.tree outdateTemp.date resultTemp.date.nexus resultTemp* temp*")

# Creating the Target Tree!
mltree <- read.tree("AssemblyMLTrees/GlobalReview.treefile")
#mltree <- read.tree("NewMLTrees/WMed.treefile")
mltree$tip.label <- gsub('[.].+', '', mltree$tip.label)
mltree$tip.label <- dates1$newname[match(mltree$tip.label, dates1$V1)]
mltree <- root(mltree, outgroup = outgroupLabel, resolve.root = T)
write.nexus(mltree, file = "GlobalREviewWithWdates.nexus")

#system('~/Applications/BEASTv1.10.4/bin/treeannotator -heights mean datedBoot.nexus datedBootSummaryTree.nexus')
#system('~/Applications/BEASTv1.10.4/bin/treeannotator -heights mean -target WMedWithWdates.nexus datedBoot.nexus WMedSummary.nexus')
system('~/Applications/BEASTv1.10.4/bin/treeannotator -heights mean -target GlobalReviewWithWdates.nexus datedBootAssemblyGlobalReview.nexus AssemblyMLTrees/GlobalReviewSummary.nexus')

cat("Rate:", quantile(rates, c(0.5, 0.025, 0.975)), "\n")
cat("tMRCA:", quantile(tMRCA, c(0.5, 0.025, 0.975)), "\n")
cat("obFunc:", quantile(obFunc, c(0.5, 0.025, 0.975)), "\n")
