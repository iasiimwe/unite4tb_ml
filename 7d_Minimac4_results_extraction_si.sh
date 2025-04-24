#!/bin/bash
# set -x  # This will print each command as it runs

for j in {1..100}; do
    for imputation_method in "michigan_afr" "michigan_eur"; do
        ref_population=$(echo "$imputation_method" | sed 's/michigan_//')

        # Initialize CSV file with headers if not already present
        csv_file="michigan/time_${ref_population}.csv"
        if [[ ! -f "$csv_file" ]]; then            
            echo "dataset,start,end" > "$csv_file"
        fi

        # Start time for the current iteration
        start=$(date +"%Y-%m-%dT%H:%M:%S.%6N")
        
        # Let us get ids
        job_id=$(awk -v j="$j" -F',' 'NR > 1 && NR == j+1 {print $2}' "michigan/job_id_${ref_population}.csv" | sed 's/"//g')

        # Check if the id was successful
        job_id_path="/pub59/iasiimwe/imputationserver/workspace/${job_id}/output/chr_1.zip"

        if [[ ! -f "$job_id_path" ]]; then
            continue
        fi

        # If imputed data exists, unzip it into a new folder
        if [[ ! -d "unzipped" ]]; then
            mkdir "unzipped"
        fi

        # Unzip the data
        /pub59/iasiimwe/IWPC/7zz x -pabc123 "$job_id_path" -o./unzipped/ -y

        # Convert to binary format
        /pub59/iasiimwe/plink1.9/plink --vcf unzipped/chr1.dose.vcf.gz --make-bed --out temp_michigan

        # Get true SNPs
        true_snps_m=$(awk 'NR > 1 {print $2}' "true_covar/true_snps_${j}_GRCh37.txt")

        # Get the SNPs that have the true covar coordinates
        to_grep=$(echo "$true_snps_m" | paste -sd'|' -)
        grep -E "$to_grep" temp_michigan.bim | awk '{print $2}' > true_snps_m.txt

        # Convert to ped format and extract only the true_snps
        /pub59/iasiimwe/plink1.9/plink --bfile temp_michigan --extract true_snps_m.txt --recode tab --out michigan/chr${j}_${job_id}_extracted 

        # End time for the current iteration
        end=$(date +"%Y-%m-%dT%H:%M:%S.%6N")
        
        # Save dataset (j), the start and end times to the CSV file
        echo "\"$j\",\"$start\",\"$end\"" >> "$csv_file"

        # Output the timing information
        echo "##################"
        echo "Method: $imputation_method"
        echo "$(awk "BEGIN {printf \"%.f\", $j}")% complete"
        echo "##################"
    done
done
