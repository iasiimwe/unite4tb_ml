#!/bin/bash
# set -x  # This will print each command as it runs

for imputation_method in "michigan_afr" "michigan_eur"; do
     ref_population=$(echo "$imputation_method" | sed 's/michigan_//')
     # Initialize a CSV file with headers 
     csv_file="michigan/info_scores_${ref_population}.csv"
     echo "dataset,CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO" > "$csv_file"

     for j in {1..100}; do          
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

        # Get true SNPs
        true_snps_m=$(awk 'NR > 1 {print $2}' "true_covar/true_snps_${j}_GRCh37.txt")
        to_filter=$(echo "$true_snps_m" | paste -sd'|' -)
        
        # Filter the data, add the dataset column, and append to the CSV file
          # .info.gz files are too large to view unzipped, so use the command gzip -dc unzipped/chr1.info.gz | less (press ‘q’ to quit)
        gzip -dc unzipped/chr1.info.gz | awk -v pattern="$to_filter" -v dataset="$j" '
           $2 ~ pattern {print dataset "," $1 "," $2 "," $3 "," $4 "," $5 "," $6 "," $7 "," $8 }' >> "$csv_file"

        # Output the timing information
        echo "##################"
        echo "Method: $imputation_method"
        echo "$(awk "BEGIN {printf \"%.f\", $j}")% complete"
        echo "##################"
    done
done
