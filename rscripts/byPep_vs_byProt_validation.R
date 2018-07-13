
library(dplyr)

byPep <- read.csv('testFiles/byPeptide_fixed.tsv', sep = '\t', stringsAsFactors = FALSE)
byProt <- read.csv('testFiles/byProtein_fixed.tsv', sep = '\t', stringsAsFactors = FALSE)

byPep_prot <- read.csv('testFiles/byPeptide_prot.tsv', sep = '\t', stringsAsFactors = FALSE)
byProt_prot <- read.csv('testFiles/byProtein_prot.tsv', sep = '\t', stringsAsFactors = FALSE)

byPep <- byPep[byPep$DatType == 'peptide',]
byProt <- byProt[byProt$DatType == 'peptide',]

matchCols <- c('ID', 'Protein', 'Description', 'Sequence', 'Charge')

noMatch <- dplyr::anti_join(byPep, byProt, by = matchCols)
