#!/bin/bash

# Usage: Obtain the identity difference between different clades of Meloidogyne

patterns=("Mg" "Mc27" "MhA1" "ME_e1834" "Mare" "Mjav" "Minc" "CACSLI")

# Function to create BLAST database
create_blast_db() {
    local fasta_file=$1
    $SING $SING_IMG makeblastdb -in "$fasta_file" -dbtype prot
}

# Function to run BLASTP
run_blastp() {
    local query_file=$1
    local db_file=$2
    local output_file=$3
    blastp -query "$query_file" -db "$db_file" -out "$output_file" -outfmt '6 qseqid sseqid pident'
}

# Function to extract percent identities based on folder type
extract_percent_identities() {
    local output_file=$1
    local folder_type=$2
    local tmp_file=$(mktemp)
    
    # Exclude same sequence comparisons
    awk '{if ($1 != $2) print $0}' "$output_file" > "$tmp_file"
    
    local percent_identities=()
    while read -r line; do
        qseqid=$(echo "$line" | cut -f1)
        sseqid=$(echo "$line" | cut -f2)
        pident=$(echo "$line" | cut -f3)

        qspecies=""
        sspecies=""

        for pattern in "${patterns[@]}"; do
            if [[ $qseqid == $pattern* ]]; then
                qspecies=$pattern
            fi
            if [[ $sseqid == $pattern* ]]; then
                sspecies=$pattern
            fi
        done

        case $folder_type in
            "clade_1")
                if [[ $qspecies == "ME_e1834" && $qspecies != $sspecies ]]; then
                    percent_identities+=("$pident")
                fi
                ;;
            "clade_1_nonment")
                if [[ $qspecies != $sspecies ]]; then
                    percent_identities+=("$pident")
                fi
                ;;
            "clade_2")
                if [[ $qspecies == "MhA1" && $qspecies != $sspecies ]]; then
                    percent_identities+=("$pident")
                fi
                ;;
            "clade_3")
                if [[ ($qspecies == "Mg" || $qspecies == "Mc27") && $qspecies != $sspecies ]]; then
                    if [[ !($qspecies == "Mg" && $sspecies == "Mc27") && !($qspecies == "Mc27" && $sspecies == "Mg") ]]; then
                        percent_identities+=("$pident")
                    fi
                fi
                ;;
        esac
    done < "$tmp_file"

    rm "$tmp_file"
    echo "${percent_identities[@]}"
}

# Function to calculate mean
calculate_mean() {
    local values=("$@")
    local sum=0
    local count=${#values[@]}
    for value in "${values[@]}"; do
        sum=$(awk -v sum="$sum" -v value="$value" 'BEGIN {print sum + value}')
    done
    echo $(awk -v sum="$sum" -v count="$count" 'BEGIN {print sum / count}')
}

# Function to calculate median
calculate_median() {
    local values=($(printf '%s\n' "${@}" | sort -n))
    local count=${#values[@]}
    if (( $count % 2 == 0 )); then
        local mid=$((count / 2))
        echo $(awk -v a="${values[$mid-1]}" -v b="${values[$mid]}" 'BEGIN {print (a + b) / 2}')
    else
        local mid=$(((count + 1) / 2))
        echo "${values[$mid-1]}"
    fi
}

# Process
process_folder() {
    local folder=$1
    local folder_type=$2
    local output_dir="output"

    mkdir -p "$output_dir"

    local output_file="$output_dir/${folder}_global_identities.txt"

    echo "Processing folder: $folder"

    mean_identities=()
    median_identities=()

    for fasta_file in "$folder"/*.fa; do
        file_name=$(basename "$fasta_file")
        db_file="$fasta_file"
        output_file_blast="$fasta_file.blastp.out"
        
        echo "Processing file: $file_name"

        create_blast_db "$fasta_file"
        run_blastp "$fasta_file" "$db_file" "$output_file_blast"
        
        percent_identities=($(extract_percent_identities "$output_file_blast" "$folder_type"))
        
        if [ ${#percent_identities[@]} -eq 0 ]; then
            mean_identity=0
            median_identity=0
        else
            mean_identity=$(calculate_mean "${percent_identities[@]}")
            median_identity=$(calculate_median "${percent_identities[@]}")
        fi

        mean_identities+=("$mean_identity")
        median_identities+=("$median_identity")

        echo "Mean identity for $file_name: $mean_identity"
        echo "Median identity for $file_name: $median_identity"
    done

    # Calculate global mean of means and mean of medians
    global_mean_of_means=$(calculate_mean "${mean_identities[@]}")
    global_mean_of_medians=$(calculate_mean "${median_identities[@]}")

    echo "Global mean of means for $folder: $global_mean_of_means" > "$output_file"
    echo "Global mean of medians for $folder: $global_mean_of_medians" >> "$output_file"
}

if [ $# -ne 1 ]; then
    echo "Usage: $0 <folder>"
    echo "Available folders: clade_1, clade_1_nonment, clade_2, clade_3"
    exit 1
fi

folder=$1

# Define folders and their types
declare -A folders
folders=(["clade_1"]="clade_1" ["clade_1_nonment"]="clade_1_nonment" ["clade_2"]="clade_2" ["clade_3"]="clade_3")

# Process the specified folder
if [ -d "$folder" ]; then
    if [ -n "${folders[$folder]}" ]; then
        process_folder "$folder" "${folders[$folder]}"
    else
        echo "Invalid folder type. Available folders: clade_1, clade_1_nonment, clade_2, clade_3"
        exit 1
    fi
else
    echo "Folder $folder does not exist"
    exit 1
fi

