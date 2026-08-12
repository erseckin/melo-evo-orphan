#!/bin/bash

# Usage: Given orphan orthogroups from both tools, identify common orthogroups

i=1
while IFS= read -r line; do
    echo "$line" | tr ',' '\n' | sort > "of_ogs/of_$i.txt"
    ((i++))
done < of.csv

j=1
while IFS= read -r line; do
    echo "$line" | tr ',' '\n' | sort > "sp_ogs/sp_$j.txt"
    ((j++))
done < sp.csv

for file in sp_ogs/*.txt; do
    while IFS= read -r line; do
        if [[ "$line" =~ [[:space:]] ]]; then
            echo "$line" | tr ' ' '\n'
        else
            echo "$line"
        fi
    done < "$file" | sort > "${file}.tmp" && mv "${file}.tmp" "$file"
done

for dir in of_ogs sp_ogs; do
    for file in "$dir"/*.txt; do
        sed -i '/^NXFT/d' "$file"  
        sed -i '/^$/d' "$file"     
    done
done

declare -A of_hashes sp_hashes

for file in of_ogs/*.txt; do
    hash=$(md5sum "$file" | awk '{print $1}')
    of_hashes[$hash]+="${file} "
done

for file in sp_ogs/*.txt; do
    hash=$(md5sum "$file" | awk '{print $1}')
    sp_hashes[$hash]+="${file} "
done

for hash in "${!of_hashes[@]}"; do
    if [[ -n "${sp_hashes[$hash]}" ]]; then
        read -ra OF_FILES <<< "${of_hashes[$hash]}"
        read -ra SP_FILES <<< "${sp_hashes[$hash]}"
        for of_file in "${OF_FILES[@]}"; do
            of_filename=$(basename "$of_file")
            for sp_file in "${SP_FILES[@]}"; do
                sp_filename=$(basename "$sp_file")
                
                cp "$of_file" "common_ogs/${of_filename}_${sp_filename}"
                break 
            done
        done
    fi
done
