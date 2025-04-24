# Load libraries
library(data.table)
library(tidyverse)

# Number of datasets
n_datasets <- 100
plink_path <- "/pub59/iasiimwe/plink1.9/plink"

# Set working directory
setwd("/pub59/iasiimwe/TB/datasets")

# Add folders if they don't exist
dirs <- c("missing", "missing/5", "missing/10", "missing/20", "missing/50")
lapply(dirs, function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE))

# Set seed for reproducibility
set.seed(7)

# Function to calculate_missingness
calculate_missingness <- function(df) {
total_elements <- prod(dim(df))  # Total number of elements in the dataframe
missing_elements <- sum(is.na(df) | df == "0 0")  # Count missing values ("NA" or "0 0")
percentage_missing <- round((missing_elements / total_elements) * 100, 2)  # Calculate percentage
return(percentage_missing)
}

# Dominant coding function
dominant_coding_fn <- function(x) {
  x <- gsub(" ", "", x) 
  y <- table(x)
  y <- names(sort(-y)) # Sort in descending order, assume mutant type is the least abundant
  if (length(y) == 3) x <- gsub(y[[3]], "1", x) # Give mutant homozygotes "1"
  if (length(y) > 1) x <- gsub(y[[2]], "1", x) # Give heterozygotes "1" 
  x <- gsub(y[[1]], "2", x) # Give homozygotes 2, so that the coding is similar to males vs females (higher 2 has less risk)
  # return(as.numeric(x))
  return(x)
}

# Missing [not] at random function
    # MAR, The total missingness is percent_missing_x, which is spread across both males and females based on their proportions.
    # Males (or mutant alelle carriers) are three times more likely to have missing data compared to females (or wild-type)
mar_mnar_fn <- function(x, percent_missing_x, missingness_factor = 3) {
  # Male and female proportions
  male_p <- prop.table(table(x))[1]
  female_p <- prop.table(table(x))[2]
  # Probability of missingness for females
  p_f <- percent_missing_x / ( missingness_factor * male_p + female_p)
  # Probability of missingness for males
  p_m <-  missingness_factor * p_f
  # Assign probabilities of missingness based on gender
  prob_missing <- ifelse(x == 1, p_m, p_f)  
  return(prob_missing)
}

# Add missing data ("MCAR", "MAR", "MNAR"; 5, 10, 20, 50%)
missingness_fn <- function(percent_missing = 0.05, mechanism = "MCAR", dat, sex) { # dat = ped_snps, sex = ped_start$V5
  # Checks
  num_snps_all <- ncol(dat)
  if (num_snps_all != 1000) stop("Error: The number of SNPs is not equal to 1000.")
  
  num_snps_true <- dat %>% select(contains("true")) %>% ncol()
  if (num_snps_true != 9) stop("Error: The number of 'true' SNPs is not equal to 9.")
  
  # False covariates (all MCAR)
  # --------------------------
  false_df <- dat %>% 
    select(!contains("true")) %>% # Get dataframe for false covariates
    mutate(ID = row_number()) %>% # Add ID column
    pivot_longer(cols = -ID, names_to = "SNP", values_to = "Genotype") # Convert to long format
  
  # Generate missingness indicator
  missing <- rbinom(nrow(false_df), size = 1, prob = percent_missing)
  
  # Add SNP values with missingness
  false_df$Genotype <- ifelse(missing == 1, "0 0", false_df$Genotype)
  
  # Convert back to wide format
  false_df <- false_df %>%
    pivot_wider(names_from = SNP, values_from = Genotype) %>%
    select(-ID)
  MCAR_false <- calculate_missingness(false_df)
  
  # True covariates
  # ---------------
  true_df <- dat %>% 
    select(contains("true")) %>% # Get dataframe for false covariates
    mutate(ID = row_number(),
           Sex = sex) %>% # Add ID column
    pivot_longer(cols = -c(ID, Sex), names_to = "SNP", values_to = "Genotype") # Convert to long format
  
  # Generate missingness indicator
  # MCAR
  if (mechanism == "MCAR") missing <- rbinom(nrow(true_df), size = 1, prob = percent_missing)
  # MAR, The total missingness is percent_missing, which is spread across both males and females based on their proportions.
  if (mechanism == "MAR") missing <- rbinom(nrow(true_df), size = 1, prob = mar_mnar_fn(true_df$Sex, percent_missing))
  # MNAR
  if (mechanism == "MNAR") {
    # Function to add missingness
    add_SNP_p_fn <- function(x) return(mar_mnar_fn(x, percent_missing))
    # Add probabilities column-wise since the SNPs have different MAFs
    missing_prob <- dat %>% 
      select(contains("true")) %>%
      mutate_all(dominant_coding_fn) %>%
      mutate_all(add_SNP_p_fn) %>%
      mutate(ID = row_number()) %>% 
      pivot_longer(cols = -ID, names_to = "SNP", values_to = "Genotype") %>%
      pull(Genotype)
    missing_prob[missing_prob >= 1] <- 0.9999999 # Set this as maximum for any values beyond 1
    missing_prob[missing_prob <= 0] <- 10^-6
    missing <- rbinom(nrow(true_df), size = 1, prob = missing_prob)
  }
  
  # Add SNP values with missingness
  true_df$Genotype <- ifelse(missing == 1, "0 0", true_df$Genotype)
  
  # Convert back to wide format
  true_df <- true_df %>%
    pivot_wider(names_from = SNP, values_from = Genotype) %>%
    select(-ID, - Sex)
  M_true <- calculate_missingness(true_df)
  
  dat <- false_df %>%
    bind_cols(true_df) %>%
    select(any_of(snps)) # Re-arrange the column names to match the original
  
  return(list(dat, MCAR_false, M_true))
}

#Results storage matrix
results_matrix <- matrix(0, nrow = n_datasets, ncol = 16,
                         dimnames = list(paste0("dataset_", 1:n_datasets), 
                                         c(paste0(c("FMCAR_", "MCAR_", "MAR_", "MNAR_"), 5),
                                           paste0(c("FMCAR_", "MCAR_", "MAR_", "MNAR_"), 10),
                                           paste0(c("FMCAR_", "MCAR_", "MAR_", "MNAR_"), 20),
                                           paste0(c("FMCAR_", "MCAR_", "MAR_", "MNAR_"), 50))))

for (i in 1:n_datasets) {
  # Get complete i'th dataset (covert to ped, tab delimited)
  system(paste0(plink_path, " --bfile 10_3/dat_10_3_", i, " --recode tab --out temp"))
  
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
  
  mechanisms <- c("MCAR", "MAR", "MNAR")
  missing_percent_values <- c(0.05, 0.1, 0.2, 0.5)
  for (j in seq_along(mechanisms)) {
    for (k in seq_along(missing_percent_values)) {
      missing_dat <- missingness_fn(percent_missing = missing_percent_values[[k]], 
                                    mechanism = mechanisms[[j]], dat = ped_snps, ped_start$V5)
      # Missing dataset ped file
      write.table(bind_cols(ped_start, missing_dat[[1]]), 
                  file = "temp.ped", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")

      # Save binary files
      system(paste0(plink_path, " --file temp --make-bed --out missing/", 
                    missing_percent_values[[k]] * 100, "/dat_", mechanisms[[j]], "_", i))

      # True covariate missingness percentages
      # Combine mechanism and missing percentage to find the column name
      col_name <- paste0(mechanisms[j], "_", missing_percent_values[k] * 100)
      # Find the column index for the respective column name
      col_index <- which(colnames(results_matrix) == col_name)
      # Assign results
      results_matrix[i, col_index] <- missing_dat[[3]]
      
      # False covariate missingness percentages
      if (mechanisms[[j]] == "MCAR") { # Will be the same for the rest, but save time by computing this once
        col_name <- paste0("F", col_name)
        col_index <- which(colnames(results_matrix) == col_name)
        results_matrix[i, col_index] <- missing_dat[[2]]
      } 
    }
  }
  message(paste0("##################\n##################\n##################\n", 
                 round(i * 100/n_datasets, 2), 
                 "% complete\n##################\n##################\n##################"))
}
write.csv(results_matrix, file = "missingness_percentages.csv")


