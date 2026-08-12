#!/bin/bash

# Usage: Eliminate orphan candidates with hits from Diamond

awk '{print $1}' diamond_lca_filtered.tsv >> non-orphans_of_diamond.txt

grep -Fvf non-orphans_of_diamond.txt Orthogroups_orphans_with_id.tsv >> Orthogroups_orphans_with_id_afterdiamond.tsv


grep -o "Mare_v4_contig_[0-9]*g[0-9]*" Orthogroups_orphans_with_id_afterdiamond.tsv > orphan_contig_list_afterdiamond.txt              # for Me.aren
grep -o "Mc27_g[0-9]*t[0-9]*_[0-9]*" Orthogroups_orphans_with_id_afterdiamond.tsv > orphan_contig_list_afterdiamond.txt                # for Me.chit
grep -o "ME_e1834_[0-9]*g[0-9]*" Orthogroups_orphans_with_id_afterdiamond.tsv >> orphan_contig_list_afterdiamond.txt                       # for Me.ente
grep -o "Mg[0-9]*.[0-9]*" Orthogroups_orphans_with_id_afterdiamond.tsv >> orphan_contig_list_afterdiamond.txt               # for Me.gram
grep -o "MhA1_Contig[0-9]*.frz[0-9]*.gene[0-9]*" Orthogroups_orphans_with_id_afterdiamond.tsv >> orphan_contig_list_afterdiamond.txt   # for Me.hapl
grep -o "Minc_v4_shac_contig_[0-9]*g[0-9]*" Orthogroups_orphans_with_id_afterdiamond.tsv >> orphan_contig_list_afterdiamond.txt        # for Me.inco
grep -o "Mjav_v4_contig_[0-9]*g[0-9]*" Orthogroups_orphans_with_id_afterdiamond.tsv >> orphan_contig_list_afterdiamond.txt             # for Me.java
grep -o "CACSLI[0-9]*.[0-9]*g[0-9]*" Orthogroups_orphans_with_id_afterdiamond.tsv >> orphan_contig_list_afterdiamond.txt               # for Me.luci
