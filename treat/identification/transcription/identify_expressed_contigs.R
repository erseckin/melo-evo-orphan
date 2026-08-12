library(tidyverse)
library(dplyr)
library(matrixStats)

# Usage: Treat RSEM FPKM matrix obtained from rnaseq to determine expressed genes
# once all expressed genes are identified, check if orphan orthogroups contain at least
# one expressed gene with grep

# M. javanica

FPKM_mjav <- read.table(file = 'matrice_rsem_FPKM_mjav.txt', sep = '\t', header = FALSE)

FPKM_mjav_p1 <- mutate(FPKM_mjav, FPKM_p1 = FPKM_mjav[,2]+1)
FPKM_mjav_p1lg <- mutate(FPKM_mjav_p1, FPKM_p1lg = log10(FPKM_mjav_p1[,3]))

FPKM_boxplot <- FPKM_mjav_p1lg[,4]
boxplot(FPKM_boxplot)
quantile(FPKM_mjav_p1lg[,4],prob=0.25,type=1)

FPKM_normalized <- FPKM_mjav_p1lg[,c(1,4)]
FPKM_normalized_filtered <- filter(FPKM_normalized,FPKM_p1lg >0.1846914)

write_tsv(FPKM_normalized_filtered, file='mjav_FPKM_expressed.tsv')

# M. arenaria

FPKM_mare <- read.table(file = 'matrice_rsem_FPKM_mare.txt', sep = '\t', header = FALSE)

FPKM_mare_p1 <- mutate(FPKM_mare, FPKM_p1 = FPKM_mare[,2]+1)
FPKM_mare_p1lg <- mutate(FPKM_mare_p1, FPKM_p1lg = log10(FPKM_mare_p1[,3]))

FPKM_boxplot_mare <- FPKM_mare_p1lg[,4]
boxplot(FPKM_boxplot_mare)
quantile(FPKM_mare_p1lg[,4],prob=0.25,type=1)

FPKM_normalized_mare <- FPKM_mare_p1lg[,c(1,4)]
FPKM_normalized_filtered_mare <- filter(FPKM_normalized_mare,FPKM_p1lg >0.1903317)

write_tsv(FPKM_normalized_filtered_mare, file='mare_FPKM_expressed.tsv')

# M. incognita

FPKM_minc <- read.table(file = 'matrice_rsem_FPKM_minc.txt', sep = '\t', header = TRUE)

FPKM_minc_with_medians <- FPKM_minc %>% mutate(Fem_median = rowMedians(as.matrix(FPKM_minc[,c(2,3,4)])), J2_median = rowMedians(as.matrix(FPKM_minc[,c(5,6,7)])),J3_median= rowMedians(as.matrix(FPKM_minc[,c(8,9,10)])), W_median= rowMedians(as.matrix(FPKM_minc[,c(11,12,13)])))
FPKM_minc_medians <- FPKM_minc_with_medians[,c(1,14,15,16,17)]

FPKM_minc_p1 <- mutate(FPKM_minc_medians, FPKM_p1_Fem = FPKM_minc_medians[,2]+1,FPKM_p1_J2 = FPKM_minc_medians[,3]+1,FPKM_p1_J3 = FPKM_minc_medians[,4]+1 ,FPKM_p1_W = FPKM_minc_medians[,5]+1)
FPKM_minc_p1lg <- mutate(FPKM_minc_p1, FPKM_p1lg_Fem = log10(FPKM_minc_p1[,6]),FPKM_p1lg_J2 = log10(FPKM_minc_p1[,7]), FPKM_p1lg_J3 = log10(FPKM_minc_p1[,8]),  FPKM_p1lg_W = log10(FPKM_minc_p1[,9]))

FPKM_normalized_minc <- FPKM_minc_p1lg[,c(1,10,11,12,13)]

boxplot(FPKM_normalized_minc[,5])
quantile(FPKM_normalized_minc[,2],prob=0.25,type=1)
quantile(FPKM_normalized_minc[,3],prob=0.25,type=1)
quantile(FPKM_normalized_minc[,4],prob=0.25,type=1)
quantile(FPKM_normalized_minc[,5],prob=0.25,type=1)

FPKM_normalized_filtered_minc <- filter(FPKM_normalized_minc,FPKM_p1lg_Fem >0 | FPKM_p1lg_J2 >0.3891661 | FPKM_p1lg_J3 >0 | FPKM_p1lg_W >0.3404441 )
write_tsv(FPKM_normalized_filtered_minc, file='minc_FPKM_expressed.tsv')

# M. graminicola

FPKM_mgra <- read.table(file = 'mgra_FPKM.tsv', sep = '\t', header = TRUE)

FPKM_mgra_with_medians <- FPKM_mgra %>% mutate(Egg_median = rowMedians(as.matrix(FPKM_mgra[,c(2,3,4)])), Fem_median = rowMedians(as.matrix(FPKM_mgra[,c(5,6,7)])),J2_median= rowMedians(as.matrix(FPKM_mgra[,c(8,9,10)])), J3_median= rowMedians(as.matrix(FPKM_mgra[,c(11,12)])))
FPKM_mgra_medians <- FPKM_mgra_with_medians[,c(1,13,14,15,16)]

FPKM_mgra_p1 <- mutate(FPKM_mgra_medians, FPKM_p1_egg = FPKM_mgra_medians[,2]+1,FPKM_p1_fem = FPKM_mgra_medians[,3]+1,FPKM_p1_J2 = FPKM_mgra_medians[,4]+1 ,FPKM_p1_J3 = FPKM_mgra_medians[,5]+1)
FPKM_mgra_p1lg <- mutate(FPKM_mgra_p1, FPKM_p1lg_egg = log10(FPKM_mgra_p1[,6]),FPKM_p1lg_fem = log10(FPKM_mgra_p1[,7]), FPKM_p1lg_J2 = log10(FPKM_mgra_p1[,8]),  FPKM_p1lg_J3 = log10(FPKM_mgra_p1[,9]))

FPKM_normalized_mgra <- FPKM_mgra_p1lg[,c(1,10,11,12,13)]

boxplot(FPKM_normalized_mgra[,2])
quantile(FPKM_normalized_mgra[,2],prob=0.25,type=1)
quantile(FPKM_normalized_mgra[,3],prob=0.25,type=1)
quantile(FPKM_normalized_mgra[,4],prob=0.25,type=1)
quantile(FPKM_normalized_mgra[,5],prob=0.25,type=1)

FPKM_normalized_filtered_mgra <- filter(FPKM_normalized_mgra,FPKM_p1lg_egg >0.7041505 | FPKM_p1lg_fem >0.469822 | FPKM_p1lg_J2 >0.4377506 | FPKM_p1lg_J3 >0.5178554 )
write_tsv(FPKM_normalized_filtered_mgra, file='mgra_FPKM_expressed.tsv')

# M. enterolobii

FPKM_ment <- read.table(file = 'ment_FPKM.tsv', sep = '\t', header = TRUE)

FPKM_ment_with_medians <- FPKM_ment %>% mutate(J2_median = rowMedians(as.matrix(FPKM_ment[,c(2,3,4)])), J3J4_median = rowMedians(as.matrix(FPKM_ment[,c(5,6,7)])))
FPKM_ment_medians <- FPKM_ment_with_medians[,c(1,8,9)]

FPKM_ment_p1 <- mutate(FPKM_ment_medians, FPKM_p1_J2 = FPKM_ment_medians[,2]+1,FPKM_p1_J3J4 = FPKM_ment_medians[,3]+1)
FPKM_ment_p1lg <- mutate(FPKM_ment_p1, FPKM_p1lg_J2 = log10(FPKM_ment_p1[,4]),FPKM_p1lg_J3J4 = log10(FPKM_ment_p1[,5]))

FPKM_normalized_ment <- FPKM_ment_p1lg[,c(1,6,7)]

boxplot(FPKM_normalized_ment[,2])
quantile(FPKM_normalized_ment[,2],prob=0.25,type=1)
quantile(FPKM_normalized_ment[,3],prob=0.25,type=1)

FPKM_normalized_filtered_ment <- filter(FPKM_normalized_ment,FPKM_p1lg_J2 >0.7846173 | FPKM_p1lg_J3J4 >0 )
write_tsv(FPKM_normalized_filtered_ment, file='ment_FPKM_expressed.tsv')

# M. luci

FPKM_mluc <- read.table(file = 'mluc_FPKM.tsv', sep = '\t', header = TRUE)

FPKM_mluc_with_medians <- FPKM_mluc %>% mutate(H_median = rowMedians(as.matrix(FPKM_mluc[,c(2,3,4)])), N_median = rowMedians(as.matrix(FPKM_mluc[,c(5,6,7)])),T_median= rowMedians(as.matrix(FPKM_mluc[,c(8,9,10)])))
FPKM_mluc_medians <- FPKM_mluc_with_medians[,c(1,11,12,13)]

FPKM_mluc_p1 <- mutate(FPKM_mluc_medians, FPKM_p1_H = FPKM_mluc_medians[,2]+1,FPKM_p1_N = FPKM_mluc_medians[,3]+1,FPKM_p1_T = FPKM_mluc_medians[,4]+1)
FPKM_mluc_p1lg <- mutate(FPKM_mluc_p1, FPKM_p1lg_H = log10(FPKM_mluc_p1[,5]),FPKM_p1lg_N = log10(FPKM_mluc_p1[,6]), FPKM_p1lg_T = log10(FPKM_mluc_p1[,7]))

FPKM_normalized_mluc <- FPKM_mluc_p1lg[,c(1,8,9,10)]

boxplot(FPKM_normalized_mluc[,4])
quantile(FPKM_normalized_mluc[,2],prob=0.25,type=1)
quantile(FPKM_normalized_mluc[,3],prob=0.25,type=1)
quantile(FPKM_normalized_mluc[,4],prob=0.25,type=1)

FPKM_normalized_filtered_mluc <- filter(FPKM_normalized_mluc,FPKM_p1lg_H >0 | FPKM_p1lg_N >0 | FPKM_p1lg_T >0 )
write_tsv(FPKM_normalized_filtered_mluc, file='mluc_FPKM_expressed.tsv')

