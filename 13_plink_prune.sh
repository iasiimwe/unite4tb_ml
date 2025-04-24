#!/bin/bash

# Define the pruning folder
PRUNING_FOLDER="pruning"

# Check if the folder exists, if not, create it
if [ ! -d "$PRUNING_FOLDER" ]; then
    echo "Creating folder: $PRUNING_FOLDER"
    mkdir -p "$PRUNING_FOLDER"
else
    echo "Folder '$PRUNING_FOLDER' already exists. Proceeding..."
fi

# Define the number of datasets per nsnp
n_datasets=100

# Define parameter combinations
params=(
    "50 5 0.1"
    "50 5 0.3"
    "50 5 0.5"
    "100 5 0.1"
    "100 5 0.3"
    "100 5 0.5"
    "100 10 0.1"
    "100 10 0.3"
    "100 10 0.5"
    "100 50 0.1"
    "100 50 0.3"
    "100 50 0.5"
    "1000 5 0.1"
    "1000 5 0.3"
    "1000 5 0.5"
    "1000 10 0.1"
    "1000 10 0.3"
    "1000 10 0.5"
    "1000 50 0.1"
    "1000 50 0.3"
    "1000 50 0.5"
    "1000 100 0.1"
    "1000 100 0.3"
    "1000 100 0.5"
)

# Define the single output file
RESULT_FILE="$PRUNING_FOLDER/pruning_results.txt"

# Create the result file with a header
echo "N_SNP J Window Step R2 Retained_SNPs" > "$RESULT_FILE"

# Loop over nsnp values (from 3 to 6)
for nsnp in {3..6}; do
    echo "Processing N_SNP = $nsnp"

    # Loop over j values (from 1 to n_datasets)
    for ((j=1; j<=n_datasets; j++)); do
        echo "Processing dataset j = $j for N_SNP = $nsnp"

        # Define input PLINK file dynamically
        input="/pub59/iasiimwe/TB/datasets/10_${nsnp}/dat_10_${nsnp}_${j}"

        # Check if the input file exists before running PLINK
        if [ ! -f "${input}.bed" ]; then
            echo "Warning: PLINK dataset ${input} not found. Skipping..."
            continue
        fi

        # Run PLINK for each parameter set
        for p in "${params[@]}"; do
            read -r window step r2 <<< "$p"

            # Define output file name
            output_prefix="$PRUNING_FOLDER/pruned_10_${nsnp}_${j}_${window}_${step}_${r2}"

            # Run PLINK pruning
            /pub59/iasiimwe/plink1.9/plink --bfile "$input" --indep-pairwise $window $step $r2 --out "$output_prefix"

            # Check if PLINK successfully generated the .prune.in file
            if [ -f "${output_prefix}.prune.in" ]; then
                retained_snps=$(wc -l < "${output_prefix}.prune.in")
            else
                retained_snps=0  # If the file doesn't exist, assume 0 SNPs retained
                echo "Warning: PLINK did not generate ${output_prefix}.prune.in"
            fi

            # Save results
            echo "10_${nsnp} $j $window $step $r2 $retained_snps" >> "$RESULT_FILE"
        done
    done
done

echo "LD pruning completed for all nsnp and j values. Results saved in $RESULT_FILE"
