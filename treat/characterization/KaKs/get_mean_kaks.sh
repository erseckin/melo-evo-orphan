#!/bin/bash

input_file="all_kaks_results.txt" # From ParaAT
output_file="mean_kaks.txt"

awk 'NR > 1 { ka += $3; ks += $4; kaks += $5; count++ }
     END {
         if (count > 0) {
             mean_ka = ka / count
             mean_ks = ks / count
             mean_kaks = kaks / count
             printf "Mean Ka: %.6f\nMean Ks: %.6f\nMean Ka/Ks: %.6f\n", mean_ka, mean_ks, mean_kaks
         } else {
             print "No data to process"
         }
     }' "$input_file" > "$output_file"

echo "Mean values calculated and saved to $output_file"

