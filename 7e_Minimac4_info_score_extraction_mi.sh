#!/bin/bash
# set -x  # This will print each command as it runs

mechanisms=("MCAR" "MAR" "MNAR")
missing_percentages=(5 10)

for imputation_method in "michigan_afr"; do
    ref_population=$(echo "$imputation_method" | sed 's/michigan_//')

    for mechanism in "${mechanisms[@]}"; do
        for missing_percentage in "${missing_percentages[@]}"; do
            # Initialize CSV file with headers if not already present
            csv_file="michigan/info_scores_${ref_population}_${mechanism}_${missing_percentage}.csv"
            csv_time_file="michigan/info_scores_time_${ref_population}_${mechanism}_${missing_percentage}.csv"

            if [[ ! -f "$csv_file" ]]; then
                echo "dataset,CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO" > "$csv_file"
            fi

            if [[ ! -f "$csv_time_file" ]]; then
                echo "dataset,start,end" > "$csv_time_file"
            fi

            for j in {1..100}; do
                # Start time for the current iteration
                start=$(date +"%Y-%m-%dT%H:%M:%S.%6N")

                # Get job ID for the dataset
                job_id=$(awk -v j="$j" -F',' 'NR > 1 && NR == j+1 {print $2}' "michigan/job_id_${mechanism}_${missing_percentage}_${ref_population}.csv" | sed 's/"//g')

                # Check if the job ID was successful
                job_id_path="/pub59/iasiimwe/imputationserver/workspace/${job_id}/output/chr_1.zip"

                if [[ ! -f "$job_id_path" ]]; then
                    continue
                fi

                # Create unzipped folder if it doesn't exist
                if [[ ! -d "unzipped" ]]; then
                    mkdir "unzipped"
                fi

                # Unzip the imputed data
                /pub59/iasiimwe/IWPC/7zz x -pabc123 "$job_id_path" -o./unzipped/ -y

                # Get true SNPs for this dataset
                true_snps_m=$(awk 'NR > 1 {print $2}' "true_covar/true_snps_${j}_GRCh37.txt")
                to_filter=$(echo "$true_snps_m" | paste -sd'|' -)

                # Filter and append dataset information to the CSV file
                gzip -dc unzipped/chr1.info.gz | awk -v pattern="$to_filter" -v dataset="$j" '
                    $2 ~ pattern {print dataset "," $1 "," $2 "," $3 "," $4 "," $5 "," $6 "," $7 "," $8 }' >> "$csv_file"

                # End time for the current iteration
                end=$(date +"%Y-%m-%dT%H:%M:%S.%6N")

                # Save dataset (j), the start, and end times to the CSV file
                echo "$j,$start,$end" >> "$csv_time_file"

                # Output progress information
                echo "##################"
                echo "Method: $imputation_method"
                echo "Mechanism: $mechanism"
                echo "Missing %: $missing_percentage"
                echo "$(awk "BEGIN {printf \"%.f\", $j}")% complete"
                echo "##################"
            done
        done
    done
done
