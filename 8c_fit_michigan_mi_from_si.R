# Load required packages
# ---------------------
library(lixoftConnectors)
initializeLixoftConnectors(software = "monolix", path = "/pub59/iasiimwe/Lixoft/MonolixSuite2023R1/", force = TRUE)
library(tidyverse)
library(data.table)
library(ggh4x)
library(mice)

select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate

# Dominant coding function 
dominant_coding_fn <- function(x, r2, y = 1) { # SNPs were excluded so 100% missingness 
  # Define an uncertainty factor based on r2 (higher r2 -> less randomness, lower r2 -> more randomness)
  uncertainty_factor <- 1 - r2  # More uncertainty for lower r2 values
  
  # Initialize new_x to store the rounded value
  new_x <- rep(NA_real_, length(x))
  
  # Incorporate missingness into the uncertainty factor 
  uncertainty_factor <- uncertainty_factor * y
  
  # Loop over each value in x
  for (i in seq_along(new_x)) {
    # Parse genotype probabilities from a string into a numeric vector
    probs <- as.numeric(strsplit(x[i], ",")[[1]]) 
    
    # Add Gaussian noise to all probability components
    noisy_probs <- probs + rnorm(3, mean = 0, sd = uncertainty_factor)
    
    # Min-max normalization: (the most negative will have zero probability)
    noisy_probs <- (noisy_probs - min(noisy_probs)) / (max(noisy_probs) - min(noisy_probs))
    
    # Normalize to ensure sum of probabilities is 1
    noisy_probs <- noisy_probs / sum(noisy_probs)
    
    # Sample a genotype (0, 1, or 2) using the adjusted probabilities
    new_x[i] <- sample(c(0, 1, 2), size = 1, prob = noisy_probs)
  }
  
  # Now apply dominant coding
  x <- as.character(new_x + 3) # We will not rely on the original wildtypes but on the frequency (consistent with earlier coding)
  y <- table(x)
  y <- names(sort(-y)) # Sort in descending order, assume mutant type is the least abundant
  if (length(y) == 3) x <- gsub(y[[3]], "1", x) # Give mutant homozygotes "1"
  if (length(y) > 1) x <- gsub(y[[2]], "1", x) # Give heterozygotes "1" 
  x <- gsub(y[[1]], "0", x) # Give homozygotes 0
  return(as.numeric(x))
} 

dominant_coding_fn2 <- function(x) {
  x <- gsub(" ", "", x) 
  x[x == "00"] <- NA
  y <- table(x)
  y <- names(sort(-y)) # Sort in descending order, assume mutant type is the least abundant
  if (length(y) == 3) x <- gsub(y[[3]], "1", x) # Give mutant homozygotes "1"
  if (length(y) > 1) x <- gsub(y[[2]], "1", x) # Give heterozygotes "1" 
  x <- gsub(y[[1]], "0", x) # Give homozygotes 0
  return(as.numeric(x))
}

# Function to generate multiple imputed datasets
generate_imputed_datasets <- function(df, r2_df, n = 5, seed = 7) {
  set.seed(seed)
  imputed_datasets <- list()
  for (i in 1:n) {
    imputed_dataset <- df
    for (col in paste0("SNP", 1:9)) {
      imputed_dataset[[col]] <- dominant_coding_fn(x = df[[col]], r2 = pull(filter(r2_df, SNP == col))) #
    }
    imputed_datasets[[i]] <- imputed_dataset
    # imputed_datasets[[i]] <- mutate_at(df, vars(!ID2), dominant_coding_fn)
  }
  return(imputed_datasets)
}

# Rubin's rule function to pool estimates from multiple imputations
rubin_fn <- function(coef_estimates, coef_ses) {
  # Number of imputations
  m <- length(coef_estimates)
  # Calculate the mean of the estimates (pooled estimate)
  qbar <- mean(coef_estimates)
  # Calculate the mean within-imputation variance
  ubar <- mean(coef_ses^2)
  # Calculate the between-imputation variance
  B <- var(coef_estimates)
  # Calculate the pooled variance
  pooled_var <- ubar + (1 + (1/m)) * B
  # Calculate the pooled standard error
  pooled_se <- sqrt(pooled_var)
  # Degrees of freedom
  df <- (m - 1) * ((ubar + (1 + 1/m) * B)/((1 + 1/m) * B))^2 # https://pubmed.ncbi.nlm.nih.gov/36346135/
  # Return the pooled estimate, relative standard error, and degrees of freedom
  return(c(
    pooled_est = qbar,
    pooled_rse = pooled_se * 100 / qbar,
    df = df
  ))
}

## Random effects (considering log distribution)
omega_squared_fn <- function(percent_CV) log((percent_CV/100)^2 + 1)

# # Without re-inflation
omega_ka <- omega_squared_fn(68.5)
omega_Cl <- omega_squared_fn(17.5)

sim_date <-"23_12_2024"
n_datasets <- 100

imputation_methods <- c("michigan_afr_ms")
mechanisms <- c("MCAR")
missing_percentages <- c(5, 10)

for (effect_scenario in c("low", "high")) {
  effect_size <- ifelse(effect_scenario == "high", 0.5, 0.15)
  
  # Add folders if they don't exist
  if (!dir.exists(paste0("missingness_", effect_scenario))) dir.create(paste0("missingness_", effect_scenario))
  
  for (k in seq_along(imputation_methods)) {
    imputation_method <- imputation_methods[[k]]
    
    for (mechanism in mechanisms) {
      for (missing_percentage in missing_percentages) {
        message(paste0("Starting\n     Method: ",imputation_method, "\n     Mechanism: ", 
                       mechanism, "\n     missing %: ", missing_percentage))
        
        # Results storage matrix
        params_for_estimation <- c("ka_theta", "ka_theta_rse", "ka_theta_df", "ka_eta", "ka_eta_rse", "ka_eta_df", 
                                   "Cl_theta", "Cl_theta_rse", "Cl_theta_df", "Cl_eta", "Cl_eta_rse", "Cl_eta_df", 
                                   "b", "b_rse", "b_df")
        results_matrix <- matrix(0, nrow = n_datasets, ncol = length(params_for_estimation),
                                 dimnames = list(paste0("dataset_", 1:n_datasets), params_for_estimation))
        time_tb <- tibble(end = vector("character", n_datasets), 
                          start = vector("character", n_datasets))
        
        for (j in 1:n_datasets) {
          # Get time of code execution start
          start <- Sys.time()
          
          # Get DV data
          dv_dat <- read_csv(paste0("sim_DV_", effect_scenario, "/sim_data", j, "_", sim_date, ".csv"), show_col_types = FALSE) %>%
            select(!SNP1:SNP9) # Remove complete SNP data
          
          # Extract reference population from imputation method
          ref_population <- gsub("michigan_|_ms", "", imputation_method)
          
          # Get job ID for this dataset
          job_id_path <- ifelse(imputation_method == "michigan_afr", 
                                paste0("michigan/job_id_", mechanism, "_", missing_percentage, "_", ref_population, ".csv"),
                                paste0("michigan/job_id_", ref_population, ".csv"))
          
          job_id <- read_csv(job_id_path, show_col_types = FALSE) %>%
            filter(dataset == j) %>%
            pull(job_id)
          
          # Check if the data path exists (skip if it doesn't)
          dat_path <- ifelse(imputation_method == "michigan_afr", 
                             paste0("michigan/chr", j, "_", job_id, "_", mechanism, "_", missing_percentage, "_extracted.GP.FORMAT"),
                             paste0("michigan/chr", j, "_", job_id, "_extracted.GP.FORMAT"))   
          
          if(!file.exists(dat_path)) next
          
          # Get the imputed data
          ped_snps <- fread(dat_path)
          colnames(ped_snps) <- gsub("_.*", "", colnames(ped_snps))
          
          # True covariates
          true_snps <- fread(paste0("true_covar/true_snps_", j, "_GRCh37.txt"), header = TRUE, sep = " ")$GRCh37_BP
          
          # Retain only the true SNPs
          ped_snps <- ped_snps[POS %in% true_snps]
          ped_snps[, SNP := paste0("SNP_", POS)]
          setcolorder(ped_snps, "SNP")  # Move SNP to the first column
          ped_snps[, c("CHROM", "POS") := NULL]  # Remove CHROM and POS columns
          
          # Retain only the first occurrence of each duplicated SNP
          # ped_snps %>% group_by(SNP) %>% slice_head(n = 1)
          ped_snps <- ped_snps[, .SD[1], by = SNP]
          
          # Reshape from wide to long (gather samples under an "ID2" column)
          ped_snps <- melt(ped_snps, id.vars = "SNP", variable.name = "ID2", value.name = "Value")
          
          # Reshape back to wide (turn SNP values into columns)
          ped_snps <- dcast(ped_snps, ID2 ~ SNP, value.var = "Value")
          
          # Retain only the first occurrence of each duplicated column and get only the true SNPs
          ped_snps <- ped_snps %>%
            select(ID2, any_of(paste0("SNP_", true_snps))) %>%
            tibble()
          colnames(ped_snps) <- c("ID2", paste0("SNP", 1:(ncol(ped_snps) - 1)))
          
          # Check if all 9 SNPs were imputed (another reason for failure - track both)
          if(ncol(ped_snps) < 10) next
          
          # Check if the data was completely imputed - it it failed, skip
          if (max(apply(ped_snps[, sapply(ped_snps, is.character)], 2, function(x) sum(is.na(x)))) > 0) next
          
          # Get R2 values
          r2_file <- ifelse(imputation_method == "michigan_afr", 
                            list(readLines(paste0("michigan/info_scores_", ref_population, "_", mechanism, "_", missing_percentage, ".csv"))),
                            list(readLines(paste0("michigan/info_scores_", ref_population, ".csv"))))[[1]]
          r2_lines <- r2_file[str_detect(r2_file, paste0("^", j, ","))]
          r2_tb <- tibble(position = gsub("SNP_", "", true_snps)) %>%
            mutate(SNP = paste0("SNP", row_number())) %>%
            left_join(tibble(position = gsub(".*,", "", str_extract(r2_lines, "^(?:[^,]+,){2}([0-9]+)")),
                             r2 = as.numeric(gsub("R2=", "", gsub(";.*", "", str_extract(r2_lines, "R2=.*")))))) %>%
            select(-position) %>%
            group_by(SNP) %>%
            slice_head(n = 1)
          
          # Obtain the multiply imputed datasets
          imputed_dat <- generate_imputed_datasets(ped_snps, r2_tb, n = missing_percentage)
          
          # Create a tibble to store results for pooling
          pooled_tb <- tibble(
            ka = rep(NA_real_, missing_percentage),
            ka_se = NA_real_,
            Cl = NA_real_,
            Cl_se = NA_real_,
            omega_ka = NA_real_,
            omega_ka_se = NA_real_,
            omega_Cl = NA_real_,
            omega_Cl_se = NA_real_,
            b = NA_real_,
            b_se = NA_real_
          )
          
          for (m in 1:missing_percentage) {
            ped_snps <- imputed_dat[[m]] %>%
              select(SNP1:SNP9, ID2) 
            
            # Find columns that have only wild-type
            snps_to_remove <- colnames(ped_snps)[apply(ped_snps[, sapply(ped_snps, is.numeric)], 2, max) == 0]
            snps_to_remove <- snps_to_remove[str_detect(snps_to_remove, "SNP")] # Not to remove ID
            
            # Join to the dv data
            dat <- ped_snps %>%
              select(!all_of(snps_to_remove)) %>%
              left_join(dv_dat) %>%
              relocate(ID:WT, .before = 1) %>%
              select(-NTIME, -ID2) 
            write.csv(dat, paste0("temp_mice_monolix_", imputation_method, "_", effect_scenario, ".csv"), row.names = FALSE)
            
            initializeLixoftConnectors(software = "monolix", path = "/pub59/iasiimwe/Lixoft/MonolixSuite2023R1/", force = TRUE)
            newProject(data = list(dataFile = paste0("temp_mice_monolix_", imputation_method, "_", effect_scenario, ".csv"),
                                   headerTypes = c("id", "time", "evid", "occ", "amount", "observation", 
                                                   "catcov", "regressor", rep("catcov", 9 - length(snps_to_remove)))),
                       modelFile = 'tb_base_vinnard.txt')
            
            # Change error model, distribution, variability and correlations
            setErrorModel(DV = "combined1")
            setIndividualParameterDistribution(ka = "logNormal", TVV = "logNormal", TVCl = "logNormal")
            setIndividualParameterVariability(id = list (ka = TRUE, TVV = FALSE, TVCl = TRUE),
                                              OCC = list (ka = FALSE, TVV = FALSE, TVCl = FALSE))
            # getIndividualParameterModel() to confirm changes
            
            # Add SNPs as covariates
            cov_formula <- paste0("setCovariateModel(TVCl = c(", 
                                  paste0(colnames(dat)[str_detect(colnames(dat), "SNP")], " = TRUE", collapse = ", "),
                                  "))")
            eval(parse(text = cov_formula))
            
            # getPopulationParameterInformation()
            # Create a dynamic list of SNP parameters
            snp_columns <- colnames(dat)[str_detect(colnames(dat), "SNP")]
            snp_parameters <- paste0("beta_TVCl_", snp_columns, "_1 = list(initialValue = effect_size, method = 'FIXED')")
            
            # Create the full parameter list dynamically
            param_list <- c(
              "ka_pop = list(initialValue = 1.31, method = 'MLE')",
              "TVV_pop = list(initialValue = 28.57, method = 'FIXED')",
              "TVCl_pop = list(initialValue = 3.52, method = 'MLE')",
              snp_parameters,  # Include dynamically generated SNP parameters
              "omega_TVCl = list(initialValue = omega_Cl, method = 'MLE')",
              "omega_ka = list(initialValue = omega_ka, method = 'MLE')",
              "a = list(initialValue = 2.41, method = 'FIXED')",
              "b = list(initialValue = 0.22, method = 'MLE')",
              "c = list(initialValue = 1, method = 'FIXED')"
            )
            
            # Combine into a single string for setPopulationParameterInformation
            param_formula <- paste0("setPopulationParameterInformation(", 
                                    paste(param_list, collapse = ", "), 
                                    ")")
            
            # Evaluate the generated formula
            eval(parse(text = param_formula))
            
            # Set tasks in scenario
            scenario <- getScenario()
            scenario$tasks = c(populationParameterEstimation = TRUE, 
                               conditionalDistributionSampling = TRUE,
                               conditionalModeEstimation = TRUE, 
                               standardErrorEstimation = TRUE, 
                               logLikelihoodEstimation = FALSE,
                               plots = FALSE)
            scenario$linearization = FALSE
            setScenario(scenario)
            runScenario()
            
            # store the estimates and s.e. in a table
            estimates_m <- tibble(param = getEstimatedStandardErrors()$stochasticApproximation$parameter,
                                  estimate = getEstimatedPopulationParameters()[names(getEstimatedPopulationParameters()) %in% getEstimatedStandardErrors()$stochasticApproximation$parameter],
                                  rses = getEstimatedStandardErrors()$stochasticApproximation$rse) %>%
              mutate(param = gsub("TV|_pop", "", param),
                     ses = estimate * rses / 100)
            
            pooled_tb$ka[m] = pull(estimates_m %>% filter(param == "ka"), estimate)
            pooled_tb$ka_se[m] = pull(estimates_m %>% filter(param == "ka"), ses)
            pooled_tb$Cl[m] = pull(estimates_m %>% filter(param == "Cl"), estimate)
            pooled_tb$Cl_se[m] = pull(estimates_m %>% filter(param == "Cl"), ses)
            pooled_tb$omega_ka[m] = pull(estimates_m %>% filter(param == "omega_ka"), estimate)
            pooled_tb$omega_ka_se[m] = pull(estimates_m %>% filter(param == "omega_ka"), ses)
            pooled_tb$omega_Cl[m] = pull(estimates_m %>% filter(param == "omega_Cl"), estimate)
            pooled_tb$omega_Cl_se[m] = pull(estimates_m %>% filter(param == "omega_Cl"), ses)
            pooled_tb$b[m] = pull(estimates_m %>% filter(param == "b"), estimate)
            pooled_tb$b_se[m] = pull(estimates_m %>% filter(param == "b"), ses)
          }
          
          # Pool the results
          estimates <- tibble(param = c("ka", "Cl", "omega_ka", "omega_Cl", "b"),
                              estimate = NA_real_,
                              rses = NA_real_,
                              df = NA_real_)
          estimates[1, "estimate"] <- rubin_fn(pooled_tb$ka, pooled_tb$ka_se)["pooled_est"]
          estimates[1, "rses"] <- rubin_fn(pooled_tb$ka, pooled_tb$ka_se)["pooled_rse"]
          estimates[1, "df"] <- rubin_fn(pooled_tb$ka, pooled_tb$ka_se)["df"]
          estimates[2, "estimate"] <- rubin_fn(pooled_tb$Cl, pooled_tb$Cl_se)["pooled_est"]
          estimates[2, "rses"] <- rubin_fn(pooled_tb$Cl, pooled_tb$Cl_se)["pooled_rse"]
          estimates[2, "df"] <- rubin_fn(pooled_tb$Cl, pooled_tb$Cl_se)["df"]
          estimates[3, "estimate"] <- rubin_fn(pooled_tb$omega_ka, pooled_tb$omega_ka_se)["pooled_est"]
          estimates[3, "rses"] <- rubin_fn(pooled_tb$omega_ka, pooled_tb$omega_ka_se)["pooled_rse"]
          estimates[3, "df"] <- rubin_fn(pooled_tb$omega_ka, pooled_tb$omega_ka_se)["df"]
          estimates[4, "estimate"] <- rubin_fn(pooled_tb$omega_Cl, pooled_tb$omega_Cl_se)["pooled_est"]
          estimates[4, "rses"] <- rubin_fn(pooled_tb$omega_Cl, pooled_tb$omega_Cl_se)["pooled_rse"]
          estimates[4, "df"] <- rubin_fn(pooled_tb$omega_Cl, pooled_tb$omega_Cl_se)["df"]
          estimates[5, "estimate"] <- rubin_fn(pooled_tb$b, pooled_tb$b_se)["pooled_est"]
          estimates[5, "rses"] <- rubin_fn(pooled_tb$b, pooled_tb$b_se)["pooled_rse"]
          estimates[5, "df"] <- rubin_fn(pooled_tb$b, pooled_tb$b_se)["df"]
          
          # Save the results
          for (parami in params_for_estimation) {
            # Determine the action based on the suffix of `parami`
            action <- case_when(
              str_detect(parami, "_theta$") ~ "theta",
              str_detect(parami, "_theta_rse") ~ "theta_rse",
              str_detect(parami, "_theta_df") ~ "theta_df",
              str_detect(parami, "_eta$") ~ "eta",
              str_detect(parami, "_eta_rse") ~ "eta_rse",
              str_detect(parami, "_eta_df") ~ "eta_df",
              str_detect(parami, "b$") ~ "b",
              str_detect(parami, "b_rse") ~ "b_rse",
              str_detect(parami, "b_df") ~ "b_df",
              TRUE ~ NA_character_ # Default case
            )
            
            if (!is.na(action)) {
              results_matrix[j, parami] <- switch(
                action,
                "theta" = pull(filter(estimates, param == gsub("_theta$", "", parami)), estimate),
                "theta_rse" = pull(filter(estimates, param == gsub("_theta_rse", "", parami)), rses),
                "theta_df" = pull(filter(estimates, param == gsub("_theta_df", "", parami)), df),
                "eta" = pull(filter(estimates, param == paste0("omega_", gsub("_eta$", "", parami))), estimate),
                "eta_rse" = pull(filter(estimates, param == paste0("omega_", gsub("_eta_rse", "", parami))), rses),
                "eta_df" = pull(filter(estimates, param == paste0("omega_", gsub("_eta_df", "", parami))), df),
                "b" = pull(filter(estimates, param == parami), estimate),
                "b_rse" = pull(filter(estimates, param == gsub("_rse", "", parami)), rses),
                "b_df" = pull(filter(estimates, param == gsub("_df", "", parami)), df)
              )
            }
          }
          
          end <- Sys.time()
          time_tb$start[j] <- paste("T", as.character(ymd_hms(start)))
          time_tb$end[j] <- paste("T", as.character(ymd_hms(end)))
          
          message(paste0("##################\n##################\n##################\n", 
                         round(j * 100/n_datasets, 2), 
                         "% complete\n     Method: ",imputation_method, "\n     Mechanism: ", 
                         mechanism, "\n     missing %: ", missing_percentage, 
                         "\n##################\n##################\n##################"))
        }
        path_to_save <- paste0("missingness_", effect_scenario, "/", imputation_method, "_", missing_percentage, "_", mechanism)
        write.csv(results_matrix, paste0(path_to_save, ".csv"), row.names = FALSE, quote = FALSE)
        write.csv(time_tb, paste0(path_to_save, "_time.csv"), row.names = FALSE)
      }
    }
  }
}


