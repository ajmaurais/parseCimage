
library(tidyr)
library(dplyr)
library(plyr)

source('~/scripts/parseCimage/rscripts/functions.R')

# read and clean byProtein data
dat_prot <- read.csv("testFiles/byProtein.tsv", sep = '\t', stringsAsFactors = FALSE)
dat_prot$dat_protType <- 'protein'
names(dat_prot)[names(dat_prot) == 'ipi'] <- 'ID'
dat_prot[is.na(dat_prot$index),]$dat_protType <- "peptide"
dat_prot$ID <- fillTheBlanks(dat_prot$ID, " ")
dat_prot$symbol <- fillTheBlanks(dat_prot$symbol, " ")
dat_prot <- dat_prot[dat_prot$dat_protType == "peptide",]
dat_prot <- dat_prot %>% dplyr::select(ID, symbol, sequence, starts_with("mr.set_"))
dat_prot <- dat_prot %>% tidyr::gather(key = Sample, value = ratio, names(dat_prot)[grepl("^mr.set_", names(dat_prot))])
dat_prot <- dat_prot[dat_prot$ratio != 0,]

#read in and clean byPeptide data
dat_pep <- read.csv('testFiles/byPeptide.tsv', sep = '\t', stringsAsFactors = FALSE)
dat_pep <- dat_pep[is.na(dat_pep$index),]
names(dat_pep)[names(dat_pep) == 'ipi'] <- 'ID'
dat_pep <- dat_pep %>% dplyr::select(ID, symbol, sequence, starts_with("mr.set_"))
dat_pep <- dat_pep %>% tidyr::gather(key = Sample, value = ratio, names(dat_pep)[grepl("^mr.set_", names(dat_pep))])
dat_pep <- dat_pep[dat_pep$ratio != 0,]

# normalize peptide ratios to 1 by Sample
dat_prot <- normalizeRatios(dat_prot)
dat_pep <- normalizeRatios(dat_pep)

sdat_prot <- dat_prot %>% dplyr::group_by(ID, symbol, Sample) %>%
  dplyr::summarise(mr = median(removeOutliers(norm_ratio)),
                   mad = meanAbsDev(removeOutliers(norm_ratio)),
                   meanr = mean(removeOutliers(norm_ratio))) %>%
  dplyr::ungroup()
sdat_prot$Sample <- substring(sdat_prot$Sample, 4)

sdat_pep <- dat_pep %>% dplyr::group_by(ID, symbol, Sample) %>%
  dplyr::summarise(mr = median(removeOutliers(norm_ratio)),
                   mad = meanAbsDev(removeOutliers(norm_ratio)),
                   meanr = mean(removeOutliers(norm_ratio))) %>%
  dplyr::ungroup()
sdat_pep$Sample <- substring(sdat_pep$Sample, 4)

prog_byPep <- read.csv('testFiles/byPeptide_fixed_long.tsv', sep = '\t', stringsAsFactors = FALSE)
prog_byProt <- read.csv('testFiles/byProtein_fixed_long.tsv', sep = '\t', stringsAsFactors = FALSE)

prog_byPep <- prog_byPep %>% dplyr::select(ID, Sample, Median_ratio, Avg_ratio, Median_avg_dev)
prog_byProt <- prog_byProt %>% dplyr::select(ID, Sample, Median_ratio, Avg_ratio, Median_avg_dev)

sdat_pep <- plyr::join(sdat_pep, prog_byPep, by = c('ID', 'Sample'))
sdat_prot <- plyr::join(sdat_prot, prog_byProt, by = c('ID', 'Sample'))

sdat_prot$mr_diff <- abs(sdat_prot$mr - sdat_prot$Median_ratio)
sdat_pep$mr_diff <- abs(sdat_pep$mr - sdat_pep$Median_ratio)

sdat_prot$meanr_diff <- abs(sdat_prot$meanr - sdat_prot$Avg_ratio)
sdat_pep$meanr_diff <- abs(sdat_pep$meanr - sdat_pep$Avg_ratio)

sdat_prot$mad_diff <- abs(sdat_prot$mad - sdat_prot$Median_avg_dev)
sdat_pep$mad_diff <- abs(sdat_pep$mad - sdat_pep$Median_avg_dev)

for(i in list(sdat_pep, sdat_prot)){
  print(paste("max mad_diff is:", max(i$mad_diff)))
  print(paste("max meanr_diff is:", max(i$meanr_diff)))
  print(paste("max mr_diff is:", max(i$mr_diff)))
}


