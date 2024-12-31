# Load libraries and get relevant functions
# -----------------------------------------
# Assumes all packages below have been installed - if not, first complete installation.
library(tidyverse) # data processing - has ggplot (data visualization), dplyr (data transformation) etc
library(data.table)
library(faux)

# Set working directory
setwd("/pub59/iasiimwe/TB/datasets")

sim_date <-"23_12_2024"

# Get complete i'th dataset (covert to ped, tab delimited)
i <- 1
system(paste0("/pub59/iasiimwe/plink1.9/plink --bfile 10_3/dat_10_3_", i, " --recode tab --out temp"))

# Process ped and map_files
ped <- tibble(fread("temp.ped", header = FALSE)) %>%
  select(ID = V1, SEX = V5)

# Simulating clinical covariates
# ------------------------------
# True covariates are the same for all datasets (1000 to 1000000) so to use the same PK datasets 
# Also sex and WT to be similar across all datasets i.e. one simulation
# For each covariate dataset, just left-join to the covariate dataset (after renaming the true covariates, to add 'True')

# # Use Vinnard et al. 2017 simulated dataset (https://pubmed.ncbi.nlm.nih.gov/29095954/) 
# Sample size to be 100 (Vinnard was 40 for occassion 1, and 24 for occassion 2 i.e. 60% returned) - only 38 in the file below
  # We will also have a sample size of 60 for occassion 2

Vinnard_dat <- read_csv("S1Table.csv", show_col_types = FALSE) %>%
  distinct(SUBID, SEX = GENDER, WT) %>% # Males were == 1, based on the percentage (55%)
  group_by(SUBID) %>%
  slice_head(n = 1) # Only baseline weight
summary(Vinnard_dat$WT)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 37.30   50.70   55.00   55.75   59.25   77.00 # Median weight different

WT_mean <- mean(Vinnard_dat$WT) # 55.74651
WT_median <-  median(Vinnard_dat$WT) # 55
WT_mean_M <- mean(Vinnard_dat %>% filter(SEX == 1) %>% pull(WT)) # 57.57368
WT_mean_F <- mean(Vinnard_dat %>% filter(SEX == 2) %>% pull(WT)) # 54.3
WT_sd_M <- sd(Vinnard_dat %>% filter(SEX == 1) %>% pull(WT)) # 11.56437
WT_sd_F <- sd(Vinnard_dat %>% filter(SEX == 2) %>% pull(WT)) # 4.600851

# Prespecified ranges
WT_range <- range(Vinnard_dat$WT) # 37.3 77.0

# initiate variables
valid_WT <- FALSE # Flag to indicate if the generated values are valid

set.seed(7)
# Start the loop
while(!valid_WT) {
  # Males
  valid <- FALSE 
  while(!valid) {
    # Simulate WT
    cont_tb_M <- rnorm_multi(n = nrow(ped %>% filter(SEX == 1)),
                             mu = WT_mean_M, # means 
                             sd = WT_sd_M, # sds 
                             varnames = "WT",
                             empirical = FALSE) %>% 
      tibble() %>% 
      mutate(ID = ped %>% filter(SEX == 1) %>% pull(ID))
    
    # Check if WT is within the prespecified range
    if (range(cont_tb_M$WT)[1] >= WT_range[1] && range(cont_tb_M$WT)[2] <= WT_range[2]) valid <- TRUE # Set valid to TRUE to exit the loop
  }
  # Females
  valid <- FALSE 
  while(!valid) {
    # Simulate WT
    cont_tb_F <- rnorm_multi(n = nrow(ped %>% filter(SEX == 2)),
                             mu = WT_mean_F, # means 
                             sd = WT_sd_F, # sds 
                             varnames = "WT",
                             empirical = FALSE) %>% 
      tibble() %>% 
      mutate(ID = ped %>% filter(SEX == 2) %>% pull(ID))
    
    # Check if WT is within the prespecified range
    if (range(cont_tb_F$WT)[1] >= WT_range[1] && range(cont_tb_F$WT)[2] <= WT_range[2]) valid <- TRUE # Set valid to TRUE to exit the loop
  }
  
  # Combine the two
  covariate_tb <- ped %>% 
    left_join(bind_rows(cont_tb_M, cont_tb_F)) %>%
    mutate_at(vars(WT), function(x) round(x, 2))
  
  median_WT <- round(median(covariate_tb$WT))

  message(paste0("Median WT is: ", median_WT))
  if (median_WT == WT_median) valid_WT <- TRUE
}
write_csv(covariate_tb, paste0("clin_cov_dataset_", sim_date, ".csv"))
message("Clinical covariate dataset identified!")
# Median WT is: 55

