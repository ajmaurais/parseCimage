
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

# normalize peptide ratios to 1 by Sample
dat_prot <- normalizeRatios(dat_prot)

# calculate summary data
sdat_prot <- dat_prot %>% dplyr::group_by(ID, symbol, Sample) %>%
  dplyr::summarise(mr = median(removeOutliers(norm_ratio)),
                   mad = meanAbsDev(removeOutliers(norm_ratio)),
                   meanr = mean(removeOutliers(norm_ratio))) %>%
  dplyr::ungroup()
sdat_prot$Sample <- substring(sdat_prot$Sample, 4)

# read py output data
# to get data in proper format, use: parseCimage -d long -o protein -f 0 <file_name>
prog_byProt <- read.csv('testFiles/byProtein_fixed_long.tsv', sep = '\t', stringsAsFactors = FALSE)
prog_byProt <- prog_byProt %>% dplyr::select(ID, Sample, Median_ratio, Avg_ratio, Median_avg_dev)

sdat_prot <- plyr::join(sdat_prot, prog_byProt, by = c('ID', 'Sample'))
sdat_prot$mr_diff <- abs(sdat_prot$mr - sdat_prot$Median_ratio)
sdat_prot$meanr_diff <- abs(sdat_prot$meanr - sdat_prot$Avg_ratio)
sdat_prot$mad_diff <- abs(sdat_prot$mad - sdat_prot$Median_avg_dev)

# print differences
print(paste("max mad_diff is:", max(sdat_prot$mad_diff)))
print(paste("max meanr_diff is:", max(sdat_prot$meanr_diff)))
print(paste("max mr_diff is:", max(sdat_prot$mr_diff)))



