# Load libraries and get relevant functions
# -----------------------------------------
# Assumes all packages below have been installed - if not, first complete installation.
library(tidyverse) # data processing - has ggplot (data visualization), dplyr (data transformation) etc
library(data.table)
library(faux)

# Set working directory
setwd("/pub59/iasiimwe/TB/datasets")

# Additive coding function
additive_coding_fn <- function(x) {
  x <- gsub(" ", "", x) 
  y <- table(x)
  y <- names(sort(-y)) # Sort in descending order, assume mutant type is the least abundant
  if (length(y) == 3) x <- gsub(y[[3]], "2", x) # Give mutant homozygotes "2"
  if (length(y) > 1) x <- gsub(y[[2]], "1", x) # Give heterozygotes "1" 
  x <- gsub(y[[1]], "0", x) # Give homozygotes 0
  # return(as.numeric(x))
  return(x)
}

sim_date <-"23_12_2024"

# Get complete i'th dataset (covert to ped, tab delimited)
i <- 1
system(paste0("/pub59/iasiimwe/plink1.9/plink --bfile 10_3/dat_10_3_", i, " --recode tab --out temp"))

# Process ped and map_files
ped <- tibble(fread("temp.ped", header = FALSE))
map <- tibble(fread("temp.map", header = FALSE))
snps <- pull(map, V2)

ped_start <- ped %>% select(V1:V6) # This selects FID, IID, father's ID, mother's ID, sex and phenotype data

snp_columns <- colnames(ped)[!colnames(ped) %in% colnames(ped_start)]
ped_snps <- ped %>%
  select(all_of(snp_columns)) 

# True covariates
true_snps <- fread(paste0("true_covar/true_snps_", i, ".txt"), header = FALSE, sep = " ")$V1
# Add true to the true SNPs
snps[snps %in% true_snps] <- paste0("true_", snps[snps %in% true_snps])
# Change the ped_snps colnames
colnames(ped_snps) <- snps # This is also the column order to use later

# Add back ID and SEX
ped_snps$ID <- ped_start$V1 

ped_snps <- ped_snps %>%
  relocate(ID, .before = 1) %>%
  mutate(SEX = ped_start$V5, .after = ID)

# Simulating clinical covariates
# ------------------------------
# True covariates are the same for all datasets (1000 to 1000000) so to use the same PK datasets 
# Also sex, FFM and WT to be similar across all datasets i.e. one simulation
# For each covariate dataset, just left-join to the covariate dataset (after renaming the true covariates, to add 'True')

# # Use Keutzer et al.simulated dataset (https://pubmed.ncbi.nlm.nih.gov/35893785/) for covariate definitions
# Keutzer_dat <- read_csv("Data S1. Simulated dataset.csv", show_col_types = FALSE) %>%
#   distinct(ID, SEX, WT, HT) %>% # based on the WT and HT means, it is likely that 0 == F, and 1 == M
#   select(-ID)
# summary(Keutzer_dat$WT)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 40.22   51.29   55.38   55.89   60.92   75.24 # Median weight different
# all sex-specific correlations for WT and HT are 1: 
# Give the Error: The correlation you specified is impossible: correlation matrix not positive definite, so use correlation = 0.9

# FFM ~ WT correlation of 0.77 from Desipramine paper
# FFM: Denti 2015 (n = 100): median 43
#      Rockwood 2016 (n = 100): IQR 38-49
#      Chirehwa 2017 (n = 61): range 28-58
# WT: Denti 2015 (n = 100): median 52, IQR 48-57
#     Mugabo 2019 (n = 51): range 38-98
# Using Wan et al. formula (https://pubmed.ncbi.nlm.nih.gov/25524443/) to calculate mean (sd) from median (IQR)
# FFM mean (sd): 43 (8)
# WT mean (sd): 52 (7)

# Green and Dufull: Predicted normal weight ratio of 17%, https://pmc.ncbi.nlm.nih.gov/articles/PMC1884581/pdf/bcp0058-0119.pdf
# round((72.7/61.8), 2) # 1.18
# 1126 males, 72.7 ± 11.9; 1121 females, 61.8 ± 9.77.

FFM_WT_means <- c(43, 52)
FFM_WT_means_M <- c(47, 56.8) # These weight means maintain the ratio of 1.18 accounting for the proportions (44 vs 45) 
FFM_WT_means_F <- c(39.8, 48.1)
FFM_WT_sds <- c(8, 7)
corr_matrix <- matrix(c(1, 0.77, 0.77, 1), nrow = 2) # Higher than 0.77

# Prespecified ranges
FFM_range <- c(28, 58)
WT_range <- c(38, 98)

# initiate variables
valid_WT <- FALSE # Flag to indicate if the generated values are valid

set.seed(7)
# Start the loop
while(!valid_WT) {
  # Males
  valid <- FALSE 
  while(!valid) {
    # Simulate FFM and WT
    cont_tb_M <- rnorm_multi(n = nrow(ped_snps %>% filter(SEX == 1)),
                             mu = FFM_WT_means_M, # means 
                             sd = FFM_WT_sds, # sds 
                             r = corr_matrix,  
                             varnames = c("FFM", "WT"),
                             empirical = FALSE) %>% 
      tibble() %>% 
      mutate(ID = ped_snps %>% filter(SEX == 1) %>% pull(ID))
    
    # Check if FFM and WT are within the prespecified ranges
    if (range(cont_tb_M$FFM)[1] >= FFM_range[1] && range(cont_tb_M$FFM)[2] <= FFM_range[2] &&
        range(cont_tb_M$WT)[1] >= WT_range[1] && range(cont_tb_M$WT)[2] <= WT_range[2]) valid <- TRUE # Set valid to TRUE to exit the loop
  }
  # Females
  valid <- FALSE 
  while(!valid) {
    # Simulate FFM and WT
    cont_tb_F <- rnorm_multi(n = nrow(ped_snps %>% filter(SEX == 2)),
                             mu = FFM_WT_means_F, # means 
                             sd = FFM_WT_sds, # sds 
                             r = corr_matrix,  
                             varnames = c("FFM", "WT"),
                             empirical = FALSE) %>% 
      tibble() %>% 
      mutate(ID = ped_snps %>% filter(SEX == 2) %>% pull(ID))
    
    # Check if FFM and WT are within the prespecified ranges
    if (range(cont_tb_F$FFM)[1] >= FFM_range[1] && range(cont_tb_F$FFM)[2] <= FFM_range[2] &&
        range(cont_tb_F$WT)[1] >= WT_range[1] && range(cont_tb_F$WT)[2] <= WT_range[2]) valid <- TRUE # Set valid to TRUE to exit the loop
  }
  
  # Combine the two
  covariate_tb <- ped_snps %>%
    select(ID, SEX) %>% 
    left_join(bind_rows(cont_tb_M, cont_tb_F)) %>%
    relocate(WT, FFM, .after = SEX) %>%
    mutate_at(vars(WT:FFM), function(x) round(x, 2))
  
  median_WT <- covariate_tb %>% group_by(SEX) %>% summarise(WT = median(WT)) %>% pull(WT)
  ratio_WT <- round(median_WT[1]/median_WT[2], 2)
  
  median_FFM <- covariate_tb %>% group_by(SEX) %>% summarise(FFM = median(FFM)) %>% pull(FFM)
  ratio_FFM <- round(median_FFM[1]/median_FFM[2], 2)
  
  median_WT <- round(median(covariate_tb$WT))
  median_FFM <- round(median(covariate_tb$FFM)) 
  
  cor_WT_FFM <- cor(covariate_tb$WT, covariate_tb$FFM) # To accept any correlation equal to or above 0.7
  
  message(paste0("WT ratio is: ", ratio_WT, "; FFM ratio is: ", ratio_FFM, "; WT is: ", median_WT, "; FFM is: ", median_FFM, "; Correlation is: ", round(cor_WT_FFM, 2)))
  if (ratio_WT >= 1.15 && ratio_WT <= 1.2 && ratio_FFM >= 1.15 && ratio_FFM <= 1.2
      && median_WT == 52 && median_FFM == 43 && cor_WT_FFM >= 0.7) valid_WT <- TRUE
}
write_csv(covariate_tb, paste0("clin_cov_dataset_", sim_date, ".csv"))
message("Clinical covariate dataset identified!")
# WT ratio is: 1.19; FFM ratio is: 1.19; WT is: 52; FFM is: 43; Correlation is: 0.73

