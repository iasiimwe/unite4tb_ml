# Load required packages (assumes a new session)
# ---------------------
library(tidyverse)
library(data.table)
library(caret)
library(glmnet) # LASSO
library(randomForest) # RF
library(lubridate) 

select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate

# Caret settings
fitControl <- trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 5)

# Dominant coding function
dominant_coding_fn <- function(x) {
  x <- gsub(" ", "", x) 
  x[x == "00"] <- NA
  y <- table(x)
  y <- names(sort(-y)) # Sort in descending order, assume mutant type is the least abundant
  if (length(y) == 3) x <- gsub(y[[3]], "1", x) # Give mutant homozygotes "1"
  if (length(y) > 1) x <- gsub(y[[2]], "1", x) # Give heterozygotes "1" 
  x <- gsub(y[[1]], "0", x) # Give homozygotes 0
  return(as.numeric(x))
}

# Get top covariates 
get_top_covar_tb_ml <- function(dataQ, mlmethod, sim_num, seedx, penalized = FALSE, select_covar = FALSE) {
  # Get importance
  if (penalized && mlmethod %in% c("lasso", "elastic_net", "ridge")) {
    if (mlmethod == "lasso") tune_alpha <- 1 else if (mlmethod == "elastic_net") tune_alpha <- 0.55 else tune_alpha <- 0
    set.seed(seedx)     
    lambdas <- train(y ~ ., data = dataQ,
                     method = "glmnet", 
                     trControl = fitControl, 
                     verbose = FALSE)$results$lambda %>%
      unique()
    trained_mod <- train(y ~ ., data = dataQ,
                         method = "glmnet", 
                         trControl = fitControl, 
                         verbose = FALSE,
                         tuneGrid = expand.grid(alpha = tune_alpha, lambda = lambdas)) 
  } else {
    set.seed(seedx)     
    trained_mod <- train(y ~ ., data = dataQ,
                         method = mlmethod, 
                         trControl = fitControl, 
                         verbose = FALSE)
  }
  # Variable importance
  top_covariates <- trained_mod %>% varImp()
  top_covar_tb <- tibble(Covariates = rownames(data.frame(top_covariates$importance)),
                         Importance = top_covariates$importance$Overall) %>%
    arrange(desc(Importance)) %>%
    filter(Importance != 0)
  top_covar_tb <- top_covar_tb[1:min(n_true_snps, nrow(top_covar_tb)), ] %>% 
    mutate(Method = mlmethod, .before = Covariates) %>% 
    mutate(Simulation = sim_num, .before = Covariates) 
  
  return(top_covar_tb)
}

# General settings
sim_date <-"23_12_2024"
n_datasets <- 100
scenarios <- c("low", "high")
n_snps <- c("10_3", "10_4", "10_5") # "10_6" possible but very computationally demanding to link PCs back to the original SNPs, which is necessary
true_covar <- paste0("SNP", 1:9)
n_true_snps <- length(true_covar)
n_round <- 1 # rounding precision (Importance)
set_seed <- 7
chosen_params <- "_1000_5_0.1.prune.in"
ml_methods <- c("glmnet", "rf") 

for (effect_scenario in scenarios) {
  # Add folders if they don't exist
  if (!dir.exists(effect_scenario)) dir.create(effect_scenario)
  folder_path <- paste0(effect_scenario, "/covariates")
  if (!dir.exists(folder_path)) dir.create(folder_path)
  
  for (n_snp in n_snps) {
    
    method_list <- rep(list(tibble()), length(ml_methods))
    time_list <- rep(list(tibble(end = vector("character", n_datasets),
                                 start = vector("character", n_datasets))), 
                     length(ml_methods))
    
    for (j in 1:n_datasets) {
      # Get DV data to get dict
      dict <- read_csv(paste0("sim_DV_", effect_scenario, "/sim_data", j, "_", sim_date, ".csv"), show_col_types = FALSE) %>%
        select(id = ID, ID2) %>%#
        distinct()
      
      # Get ETAs (to use as outcome)
      etas <- read_csv(paste0(effect_scenario, "/base_model/sim_", j, "_", "etas.csv"), show_col_types = FALSE) %>%
        filter(OCC == 1) %>%
        left_join(dict) %>%
        select(ID2, y = eta_TVCl) 

      # True covariates - we need to put them as the first 9 covariates for later tracking
      true_snps <- fread(paste0("true_covar/true_snps_", j, ".txt"), header = FALSE, sep = " ")$V1

      # Dictionary
      dict_snps <- tibble(old = true_snps) %>%
        mutate(new = paste0("SNP", row_number()))
      
      # Get pruned SNP data
      dat_patho <- paste0(n_snp, "/dat_",n_snp, "_", j)
      temp_path <- paste0("pca_o_temp_", n_snp)
      system(paste0("/pub59/iasiimwe/plink1.9/plink --bfile ", dat_patho, " --recode tab --out ", temp_path))
      ped <- tibble(fread(paste0(temp_path, ".ped"), header = FALSE))
      map <- tibble(fread(paste0(temp_path, ".map"), header = FALSE))
      snps <- pull(map, V2)
      ped_start <- ped %>% select(V1:V6) # This selects FID, IID, father's ID, mother's ID, sex and phenotype data
      snp_columns <- colnames(ped)[!colnames(ped) %in% colnames(ped_start)]
      ped_snps <- ped %>%
        select(all_of(snp_columns)) 
      colnames(ped_snps) <- snps 
      
      ped_snps <- ped_snps %>%
        select(any_of(true_snps), everything()) %>% # This allows the order to be the same as in 'true_snps'
        data.table()
      ped_snps[, (colnames(ped_snps)) := lapply(.SD, dominant_coding_fn), .SDcols = colnames(ped_snps)]
      ped_snps$ID2 <- ped_start$V1
      
      # True covariates - we need to put them as the first 9 covariates for later tracking
      first_n_snps <- tibble(old = colnames(ped_snps)[1:n_true_snps]) %>% # The true SNPs should be among the first n snps
        left_join(dict_snps) %>%
        mutate(new = ifelse(is.na(new), old, new)) %>%
        pull(new)
      colnames(ped_snps)[1:n_true_snps] <- first_n_snps
      
      # Add the outcome etas
      ped_snps <- ped_snps %>%
        left_join(etas) %>%
        relocate(y, .before = 1) %>%
        select(-ID2)

      for (k in seq_along(ml_methods)) {
        message(paste0("Starting Method: ", ml_methods[[k]]))
        
        # Get time of code execution start
        start <- Sys.time()
        
        study <- ped_snps
        
        # Perform PCA
        pca_o_result <- prcomp(select(study, -y), center = TRUE, scale. = TRUE)
        analysis_tb <- tibble(as.data.frame(pca_o_result$x)) %>%
          mutate(y = study$y, .before = 1)
        
        # Find the top-m contributing SNPs
        top_SNPs <- tibble(PC = colnames(pca_o_result$rotation),
                           top_SNPs = "")
        for (PC in seq_along(colnames(pca_o_result$rotation))) {
          pc_loadings <- rev(sort(abs(pca_o_result$rotation[, PC])))[1:n_true_snps]
          top_SNPs$top_SNPs[PC] <- paste0(names(pc_loadings), collapse = ", ")
        }
        
        # Run model (with PCAs)
        set.seed(set_seed)     
        trained_mod <- train(y ~ ., data = analysis_tb,
                             method = ml_methods[[k]], 
                             trControl = fitControl, 
                             verbose = FALSE)
        
        # PCA importance 
        # We get the top 9 PCs in case each SNP has one true SNP
        # For each PCA, we get top 9 SNPs in case only one PCA is important
        # This gives us a total of 81 SNPs that we regress again
        top_pcs <- trained_mod %>% varImp()
        top_pcs_tb <- tibble(PC = rownames(data.frame(top_pcs$importance)),
                             Importance = top_pcs$importance$Overall) %>%
          arrange(desc(Importance)) %>%
          mutate(Importance = round(Importance, 3)) %>% # Eliminates importance very close to zero
          filter(Importance != 0)
        top_pcs_tb <- top_pcs_tb[1:min(n_true_snps, nrow(top_pcs_tb)), ] %>%
          left_join(top_SNPs) %>%
          select(-Importance)
        
        # Split the comma-separated SNPs into a vector
        SNPs_to_regress <- top_pcs_tb %>%
          separate_rows(top_SNPs, sep = ", ") %>%
          pull(top_SNPs) %>%
          unique()

        if (length(SNPs_to_regress) > n_true_snps) {
          # Save SNPs to extract
          tibble(new = SNPs_to_regress) %>% 
            left_join(dict_snps) %>%
            mutate(old = ifelse(is.na(old), new, old)) %>%
            select(old) %>%
            write.table(file = paste0("pca_o_SNPs_to_extract_", n_snp, ".txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
          
          # Get new dataset
          temp_path <- paste0("pca_o_new_temp_", n_snp)
          system(paste0("/pub59/iasiimwe/plink1.9/plink --bfile ", dat_patho,
                        " --extract pca_o_SNPs_to_extract_", n_snp, ".txt --recode tab --out ",
                        temp_path))
          
          # Process ped and map_files
          ped_new <- tibble(fread(paste0(temp_path, ".ped"), header = FALSE))
          ped_new <- ped_new %>% 
            select(!V1:V6) %>%
            mutate_all(dominant_coding_fn) %>%
            mutate(ID2 = ped_new$V1, .before = 1) %>%
            left_join(etas) %>%
            relocate(y, .before = 1) %>%
            select(-ID2)
          snp_names <- tibble(fread(paste0(temp_path, ".map"), header = FALSE)) %>%
            select(old = V2) %>%
            left_join(dict_snps) %>%
            mutate(new = ifelse(is.na(new), old, new)) %>%
            pull(new)
          colnames(ped_new) <- c("y", snp_names)
          
          top_covariate_tibble <- get_top_covar_tb_ml(ped_new, ml_methods[[k]], j, set_seed)
          
        } else {
          top_covariate_tibble <- tibble(old = SNPs_to_regress) %>%
            left_join(dict_snps) %>%
            mutate(new = if_else(is.na(new), old, new)) %>%
            select(Covariates = new) %>%
            mutate(Method = ml_methods[[k]], Simulation = j, .before = Covariates)
        }
        
        method_list[[k]] <- bind_rows(method_list[[k]], top_covariate_tibble)
        end <- Sys.time()
        time_list[[k]]$start[j] <- paste("T", as.character(ymd_hms(start)))
        time_list[[k]]$end[j] <- paste("T", as.character(ymd_hms(end)))
      }
      
      message(paste0("Effect Size: ", effect_scenario, 
                     "; NSNPs: ", n_snp, 
                     "; ", round(j*100/n_datasets, 3), "% complete!"))
    }
    for (k in seq_along(ml_methods)) {
      write.csv(time_list[[k]], paste0(folder_path, "/pca_o_", n_snp, "_", ml_methods[[k]], "_time.csv"), row.names = FALSE)
      write.csv(method_list[[k]], paste0(folder_path, "/pca_o_", n_snp, "_", ml_methods[[k]], "_sim_results.csv"), row.names = FALSE)
    } 
  }
}

