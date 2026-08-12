library(GeneBridge)
library(ape)
library(tibble)

# Usage: run bridge on 8 melo species as reference

# Input

species_tree <- read.tree(file = "/path/to/species_tree.nwk")
orthogroups <- read.csv(file = "/path/to/protein_ssp_og.tsv", sep='\t')

## M.graminicola

gbr_mgra <- newBridge(ogdata=orthogroups, phyloTree=species_tree, refsp='Megram')
gbr_mgra <- runBridge(gbr_mgra, threshold = 0.3)
gbr_mgra <- runPermutation(gbr_mgra, nPermutations=1000)
res_mgra <- getBridge(gbr_mgra, what = "results")

res_mgra <- tibble::rownames_to_column(res_mgra, "OG")
write.table(res_mgra, file='bridge_mgra.tsv', quote = FALSE, row.names = FALSE, sep='\t')

## M. chitwoodi

gbr_mchit <- newBridge(ogdata=orthogroups, phyloTree=species_tree, refsp='Mechit')
gbr_mchit <- runBridge(gbr_mchit, threshold = 0.3)
gbr_mchit <- runPermutation(gbr_mchit, nPermutations=1000)
res_mchit <- getBridge(gbr_mchit, what = "results")

res_mchit <- tibble::rownames_to_column(res_mchit, "OG")
write.table(res_mchit, file='bridge_mchit.tsv', quote = FALSE, row.names = FALSE, sep='\t')

## M. hapla

gbr_mhap <- newBridge(ogdata=orthogroups, phyloTree=species_tree, refsp='Mehapl')
gbr_mhap <- runBridge(gbr_mhap, threshold = 0.3)
gbr_mhap <- runPermutation(gbr_mhap, nPermutations=1000)
res_mhap <- getBridge(gbr_mhap, what = "results")

res_mhap <- tibble::rownames_to_column(res_mhap, "OG")
write.table(res_mhap, file='bridge_mhap.tsv', quote = FALSE, row.names = FALSE, sep='\t')

## M. enterolobii

gbr_ment <- newBridge(ogdata=orthogroups, phyloTree=species_tree, refsp='Meente')
gbr_ment <- runBridge(gbr_ment, threshold = 0.3)
gbr_ment <- runPermutation(gbr_ment, nPermutations=1000)
res_ment <- getBridge(gbr_ment, what = "results")

res_ment <- tibble::rownames_to_column(res_ment, "OG")
write.table(res_ment, file='bridge_ment.tsv', quote = FALSE, row.names = FALSE, sep='\t')

## M. arenaria

gbr_mare <- newBridge(ogdata=orthogroups, phyloTree=species_tree, refsp='Mearen')
gbr_mare <- runBridge(gbr_mare, threshold = 0.3)
gbr_mare <- runPermutation(gbr_mare, nPermutations=1000)
res_mare <- getBridge(gbr_mare, what = "results")

res_mare <- tibble::rownames_to_column(res_mare, "OG")
write.table(res_mare, file='bridge_mare.tsv', quote = FALSE, row.names = FALSE, sep='\t')

## M. javanica

gbr_mjav <- newBridge(ogdata=orthogroups, phyloTree=species_tree, refsp='Mejava')
gbr_mjav <- runBridge(gbr_mjav, threshold = 0.3)
gbr_mjav <- runPermutation(gbr_mjav, nPermutations=1000)
res_mjav <- getBridge(gbr_mjav, what = "results")

res_mjav <- tibble::rownames_to_column(res_mjav, "OG")
write.table(res_mjav, file='bridge_mjav.tsv', quote = FALSE, row.names = FALSE, sep='\t')

## M. incognita

gbr_minc <- newBridge(ogdata=orthogroups, phyloTree=species_tree, refsp='Meinco')
gbr_minc <- runBridge(gbr_minc, threshold = 0.3)
gbr_minc <- runPermutation(gbr_minc, nPermutations=1000)
res_minc <- getBridge(gbr_minc, what = "results")

res_minc <- tibble::rownames_to_column(res_minc, "OG")
write.table(res_minc, file='bridge_minc.tsv', quote = FALSE, row.names = FALSE, sep='\t')

## M.luci

gbr_mluc <- newBridge(ogdata=orthogroups, phyloTree=species_tree, refsp='Meluci')
gbr_mluc <- runBridge(gbr_mluc, threshold = 0.3)
gbr_mluc <- runPermutation(gbr_mluc, nPermutations=1000)
res_mluc <- getBridge(gbr_mluc, what = "results")

res_mluc <- tibble::rownames_to_column(res_mluc, "OG")
write.table(res_mluc, file='bridge_mluc.tsv', quote = FALSE, row.names = FALSE, sep='\t')
