#!/bin/bash

n_datasets=100
mechanisms=("MCAR" "MAR" "MNAR")
missing_percentages=(5 10)

for mechanism in "${mechanisms[@]}"; do
    for missing_percentage in "${missing_percentages[@]}"; do
        for j in $(seq 1 $n_datasets); do
            for imputation_method in "michigan_afr"; do
                # Extract reference population from imputation method
                ref_population=$(echo "$imputation_method" | sed 's/michigan_//')

                # Define CSV file with headers if not already present
                csv_file="michigan/time_${mechanism}_${missing_percentage}_${ref_population}.csv"

                if [[ ! -f "$csv_file" ]]; then
                    echo "dataset,start,end" > "$csv_file"
                fi

                # Start time for the current iteration
                start=$(date +"%Y-%m-%dT%H:%M:%S.%6N")

                # Get job ID for this dataset
                job_id=$(awk -v j="$j" -F',' 'NR > 1 && NR == j+1 {print $2}' "michigan/job_id_${mechanism}_${missing_percentage}_${ref_population}.csv" | sed 's/"//g')

                # Check if the job ID path exists
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

                # Get true SNPs for the dataset
                true_snps_m=$(awk 'NR > 1 {print $2}' "true_covar/true_snps_${j}_GRCh37.txt")

                # Extract SNPs with the correct coordinates
                to_grep=$(echo "$true_snps_m" | paste -sd'|' -)
                grep -E "$to_grep" temp_michigan.bim | awk '{print $2}' > true_snps_m.txt

                # Convert to PED format, extracting only true SNPs
                # /pub59/iasiimwe/plink1.9/plink --bfile temp_michigan --extract true_snps_m.txt --recode tab --out michigan/chr${j}_${job_id}_${mechanism}_${missing_percentage}_extracted

                # Obtain genotype probabilities
                /pub59/iasiimwe/vcftools_0.1.13/bin/vcftools --gzvcf unzipped/chr1.dose.vcf.gz --out michigan/chr${j}_${job_id}_${mechanism}_${missing_percentage}_extracted --snps true_snps_m.txt --extract-FORMAT-info GP

                # End time for the current iteration
                end=$(date +"%Y-%m-%dT%H:%M:%S.%6N")

                # Save dataset, start, and end times to the CSV file
                echo "\"$j\",\"$start\",\"$end\"" >> "$csv_file"

                # Output the timing information
                echo "##################"
                echo "Method: $imputation_method | Mechanism: $mechanism | Missing: ${missing_percentage}%"
                echo "$(awk "BEGIN {printf \"%.f\", ($j/$n_datasets)*100}")% complete"
                echo "##################"
            done
        done
    done
done
