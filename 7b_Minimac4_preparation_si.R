# Load libraries
library(data.table)
library(tidyverse) 

# Add folders if they don't exist
if (!dir.exists("michigan")) dir.create("michigan")

# Relevant function
not_scientific <- function(x) trimws(format(x, scientific = FALSE))

n_datasets <- 100
# Times storage matrices
time_tb <- tibble(end = vector("character", n_datasets), 
                  start = vector("character", n_datasets))
for (k in 1:n_datasets) {
  # Get time of code execution start
  start <- Sys.time()
  
  # Paths
  true_snps <- paste0("true_covar/true_snps_", k, ".txt")
  geno_path <- paste0("10_3/dat_10_3_", k)
  
  # Convert true_snps to GRCh37 (will need these names later)
  tibble(fread(true_snps, header = FALSE, sep = "\t")) %>%
    mutate(X1 = "chr1", X2 = as.numeric(gsub("1:|\\[.*", "", V1)), X3 = X2 + 1) %>%
    mutate_at(vars(X2,X3), not_scientific) %>%
    select(X1, X2, X3, V1)  %>%
    write.table(file = "true_snps_ucsc.bim", row.names = FALSE, col.names = FALSE, quote = FALSE) 
  system("/pub59/iasiimwe/IWPC/liftOver true_snps_ucsc.bim /pub59/iasiimwe/imputationserver/hg38ToHg19.over.chain.gz true_snps_lifted.bim true_snps_unlifted.bim")
  tibble(fread("true_snps_lifted.bim", header = FALSE, sep = "\t")) %>%
    select(SNP = V4, GRCh37_BP = V2) %>%
    write.table(file = paste0("true_covar/true_snps_", k, "_GRCh37.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  # QC process, Step 2: Convert to ped/map while excluding the true_snps from the genetic data
  system(paste0("/pub59/iasiimwe/plink1.9/plink --bfile ", geno_path, " --exclude ", true_snps, " --recode --out temp_michigan"))
  
  # QC process, Step 2b: Use liftover to change to GRCh37 (https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver for chain.gz file, QC tools work with GRCh37)                         
  read_delim("temp_michigan.map", show_col_types = FALSE, col_names = FALSE) %>%
    select(X1, X4, X2) %>%
    mutate(X1 = paste0("chr", X1), X3 = X4 + 1) %>%
    mutate_at(vars(X4,X3), not_scientific) %>%
    relocate(X2, .after = X3) %>%
    write.table(file = "temp_michigan_ucsc.bim", row.names = FALSE, col.names = FALSE, quote = FALSE)     
  system("/pub59/iasiimwe/IWPC/liftOver temp_michigan_ucsc.bim /pub59/iasiimwe/imputationserver/hg38ToHg19.over.chain.gz temp_michigan_lifted.bim temp_michigan_unlifted.bim")
  
  # Change the positions based on the lifted over coordinates
  bim_new <- read_delim("temp_michigan_lifted.bim", show_col_types = FALSE, col_names = FALSE) %>%
    select(SNP = X4, BP = X2)

  read_delim("temp_michigan.map", show_col_types = FALSE, col_names = FALSE) %>%
    rename(SNP = X2) %>%
    left_join(bim_new) %>%
    mutate(BP = not_scientific(if_else(is.na(BP), X4, BP))) %>%
    select(-X4) %>%
    write.table(file = "temp_michigan.map", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  
  # Convert to binary file and exclude SNPs not lifted
  system("/pub59/iasiimwe/plink1.9/plink --file temp_michigan --make-bed --out temp_michigan")
  
  # QC process, Step 3: Create a frequency file
  system("/pub59/iasiimwe/plink1.9/plink --freq --bfile temp_michigan --out temp_geno_freq")
  
  # QC process, Step 4: Use perl script
  system("/bin/perl HRC-1000G-check-bim.pl -b temp_michigan.bim -f temp_geno_freq.frq -r 1000GP_Phase3_combined.legend.gz -g -p AFR")
  
  # QC process, Step 5a: Edit scripts produced in the Step-above
  system("sed -i 's/\\bplink\\b/\\/pub59\\/iasiimwe\\/plink1.9\\/plink/g' Run-plink.sh") # Replace plink with /pub59/iasiimwe/plink1.9/plink
  system("sed -i '/--chr [^1]/d' Run-plink.sh") # Remove all lines with --chr x where x is not 1
  system("sed -i '/--chr 1[0-9]/d' Run-plink.sh") # Remove all lines with '--chr 1x' where x is a digit from 0 to 9
  
  # QC process, Step 5b: Run scripts produced in the Step-above - Perform SNP QC and convert to .vcf format
  system("/bin/chmod 755 Run-plink.sh") # Make it executable
  system("./Run-plink.sh") # Run it
  
  # QC process, Step 6: Change to vcf and bgzip (also save)
  system(paste0("/pub59/iasiimwe/plink1.9/plink --bfile temp_michigan-updated-chr1 --recode vcf --out michigan/chr", k))
  system(paste0("/pub59/iasiimwe/htslib/htslib-1.3.2/bgzip michigan/chr", k, ".vcf"))
  
  end <- Sys.time()
  
  time_tb$start[k] <- paste("T", as.character(ymd_hms(start)))
  time_tb$end[k] <- paste("T", as.character(ymd_hms(end)))
  
  message(paste0("##################\n##################\n##################\n", 
                 round(k * 100/n_datasets, 2), "% complete!", 
                 "\n##################\n##################\n##################"))
}
write.csv(time_tb, "michigan/time.csv", row.names = FALSE)

