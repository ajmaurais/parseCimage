
library(tidyr)
library(dplyr)
library(plyr)

source('~/scripts/parseCimage/rscripts/functions.R')

#read in and clean byPeptide data
dat_pep <- read.csv('testFiles/combined_dta.txt', sep = '\t', stringsAsFactors = FALSE)
dat_pep <- dat_pep[is.na(dat_pep$index),]
names(dat_pep)[names(dat_pep) == 'ipi'] <- 'ID'
dat_pep <- dat_pep %>% dplyr::select(ID, symbol, sequence, starts_with("mr.set_"))
dat_pep <- dat_pep %>% tidyr::gather(key = Sample, value = ratio, names(dat_pep)[grepl("^mr.set_", names(dat_pep))])
dat_pep <- dat_pep[dat_pep$ratio != 0,]

# normalize peptide ratios to 1 by Sample
dat_pep <- normalizeRatios(dat_pep)

# calcluate summary data
sdat_pep <- dat_pep %>% dplyr::group_by(ID, symbol, Sample) %>%
  dplyr::summarise(mr = median(removeOutliers(norm_ratio)),
                   mad = meanAbsDev(removeOutliers(norm_ratio)),
                   meanr = mean(removeOutliers(norm_ratio))) %>%
  dplyr::ungroup()
sdat_pep$Sample <- substring(sdat_pep$Sample, 4)

# read py program output
# to generate prog output data in proper format use: parseCimage -d long -o protein -f 0 -i byPeptide <file_name>
prog_byPep <- read.csv('testFiles/combined_dta_fixed_long.tsv', sep = '\t', stringsAsFactors = FALSE)
prog_byPep <- prog_byPep %>% dplyr::select(ID, Sample, Median_ratio, Avg_ratio, Median_avg_dev)

# combine two datasets and calculate differences
sdat_pep <- plyr::join(sdat_pep, prog_byPep, by = c('ID', 'Sample'))
sdat_pep$mr_diff <- abs(sdat_pep$mr - sdat_pep$Median_ratio)
sdat_pep$meanr_diff <- abs(sdat_pep$meanr - sdat_pep$Avg_ratio)
sdat_pep$mad_diff <- abs(sdat_pep$mad - sdat_pep$Median_avg_dev)

# print differences
print(paste("max mad_diff is:", max(sdat_pep$mad_diff)))
print(paste("max meanr_diff is:", max(sdat_pep$meanr_diff)))
print(paste("max mr_diff is:", max(sdat_pep$mr_diff)))


