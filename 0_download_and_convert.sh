#!/bin/bash

# Base URL for 1000 Genomes data
BASE_URL="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL"

# Create a folder to store the data
DATA_FOLDER="Gdata"
# mkdir -p ${DATA_FOLDER}

# Loop through chromosomes (only 1 chromosome but the code below allows for multiple chomosome downloads)
for CHR in {1..1}; do
    # Define file names
    VCF_FILE="ALL.chr${CHR}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
    TBI_FILE="${VCF_FILE}.tbi"
    OUTPUT_PREFIX="${DATA_FOLDER}/chr${CHR}"

    echo "Downloading VCF and index files for chromosome ${CHR}..."
    
    # Download the VCF file and its index
    wget -c "${BASE_URL}/${VCF_FILE}" -O ${DATA_FOLDER}/${VCF_FILE}
    wget -c "${BASE_URL}/${TBI_FILE}" -O ${DATA_FOLDER}/${TBI_FILE}
    
    echo "Converting chromosome ${CHR} to PLINK format..."
    
    # Convert to PLINK binary format
    /pub59/iasiimwe/plink2 --vcf ${DATA_FOLDER}/${VCF_FILE} \
          --double-id \
          --max-alleles 2 \
          --keep luhya100.txt \
          --geno 0.05 \
          --maf 0.01 \
          --hwe 0.000001 \
          --make-bed \
          --out ${OUTPUT_PREFIX}
    
    echo "Chromosome ${CHR} processing complete!"
done

echo "All chromosomes processed! Files are stored in the ${DATA_FOLDER} folder."
