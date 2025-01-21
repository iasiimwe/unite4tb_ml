# Load required packages
# ---------------------
library(lixoftConnectors)
initializeLixoftConnectors(software = "monolix", path = "/pub59/iasiimwe/Lixoft/MonolixSuite2024R1/", force = TRUE)
library(tidyverse)
library(data.table)
library(ggh4x)
library(mice)

select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate

# Add folders if they don't exist
if (!dir.exists("missingness")) dir.create("missingness")

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

# Clean up function (in case we have any values below 0 or above 2)
dominant_coding_fn <- function(x) {
  x <- round(x)
  x[x < 0] <- 0
  x[x > 1] <- 1
  return(x)
}

# Rubin's rule function to pool estimates from multiple imputations
rubin_fn <- function(coef_estimates, coef_ses) {
  # Calculate the mean of the estimates (pooled estimate)
  qbar <- mean(coef_estimates)
  # Calculate the mean within-imputation variance
  ubar <- mean(coef_ses^2)
  # Calculate the between-imputation variance
  B <- var(coef_estimates)
  # Calculate the pooled variance
  pooled_var <- ubar + (1 + (1/length(coef_estimates))) * B
  # Calculate the pooled standard error
  pooled_se <- sqrt(pooled_var)
  # Return the pooled estimate, standard error, t-statistic, degrees of freedom, and p-value
  return(c(
    pooled_est = qbar,
    pooled_rse = pooled_se * 100 / qbar
  ))
}

## Random effects (considering log distribution)
omega_squared_fn <- function(percent_CV) log((percent_CV/100)^2 + 1)

# # Without re-inflation
omega_ka <- omega_squared_fn(68.5)
omega_Cl <- omega_squared_fn(17.5)

sim_date <-"23_12_2024"
n_datasets <- 100
effect_size <- 0.5

imputation_methods <- c("pmm", "midastouch", "sample", "cart", "rf", "mean", "norm", "ri")
mechanisms <- c("MCAR", "MAR", "MNAR", "MNAR_post")
missing_percentages <- c(5, 10, 20, 50)

for (k in seq_along(imputation_methods)) {
  imputation_method <- imputation_methods[[k]]
  
  for (mechanism in mechanisms) {
    for (missing_percentage in missing_percentages) {
      message(paste0("Starting\n     Method: ",imputation_method, "\n     Mechanism: ", 
                     mechanism, "\n     missing %: ", missing_percentage))
      
      # Results storage matrix
      params_for_estimation <- c("ka_theta", "ka_theta_rse", "ka_eta", "ka_eta_rse", 
                                 "Cl_theta", "Cl_theta_rse", "Cl_eta", "Cl_eta_rse", "b", "b_rse")
      results_matrix <- matrix(0, nrow = n_datasets, ncol = length(params_for_estimation),
                               dimnames = list(paste0("dataset_", 1:n_datasets), params_for_estimation))
      time_tb <- tibble(end = vector("character", n_datasets), 
                        start = vector("character", n_datasets))
      
      for (j in 1:n_datasets) {
        # Get time of code execution start
        start <- Sys.time()
        
        # Get DV data
        dv_dat <- read_csv(paste0("sim_DV/sim_data", j, "_", sim_date, ".csv"), show_col_types = FALSE) %>%
          select(!SNP1:SNP9) # Remove complete SNP data
        
        # Path to imputed dataset
        mice_path <- paste0("mice/", imputation_method, "_", missing_percentage, "_", mechanism, "_", j, ".rds")
        if(!file.exists(mice_path)) next
        
        # Read imputed dataset
        imputed_dat <- read_rds(mice_path)
        
        # Check if the data was completely imputed - it it failed, skip
        dat_check <- tibble(complete(imputed_dat, 1, include = FALSE)) %>%
          select(SNP1:SNP9)
        if (max(apply(dat_check[, sapply(dat_check, is.numeric)], 2, function(x) sum(is.na(x)))) > 0) next
        
        # Create a tibble to store pooled results
        pooled_tb <- tibble(
          ka = rep(NA_real_, imputed_dat$m),
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

        for (m in 1:imputed_dat$m) {
          ped_snps <- tibble(complete(imputed_dat, m, include = FALSE)) %>%
            select(SNP1:SNP9) %>%
            mutate_all(dominant_coding_fn) %>%
            mutate(ID2 = unique(dv_dat$ID2))
          
          # Find columns that have only wild-type
          snps_to_remove <- colnames(ped_snps)[apply(ped_snps[, sapply(ped_snps, is.numeric)], 2, max) == 0]
          snps_to_remove <- snps_to_remove[str_detect(snps_to_remove, "SNP")] # Not to remove ID
          
          # Join to the dv data
          dat <- ped_snps %>%
            select(!all_of(snps_to_remove)) %>%
            left_join(dv_dat) %>%
            relocate(ID:WT, .before = 1) %>%
            select(-NTIME, -ID2) 
          write.csv(dat, paste0("temp_mice_monolix_", imputation_method, ".csv"), row.names = FALSE)
          
          initializeLixoftConnectors(software = "monolix", path = "/pub59/iasiimwe/Lixoft/MonolixSuite2024R1/", force = TRUE)
          newProject(data = list(dataFile = paste0("temp_mice_monolix_", imputation_method, ".csv"),
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
                            rses = NA_real_)
        estimates[1, "estimate"] <- rubin_fn(pooled_tb$ka, pooled_tb$ka_se)["pooled_est"]
        estimates[1, "rses"] <- rubin_fn(pooled_tb$ka, pooled_tb$ka_se)["pooled_rse"]
        estimates[2, "estimate"] <- rubin_fn(pooled_tb$Cl, pooled_tb$Cl_se)["pooled_est"]
        estimates[2, "rses"] <- rubin_fn(pooled_tb$Cl, pooled_tb$Cl_se)["pooled_rse"]
        estimates[3, "estimate"] <- rubin_fn(pooled_tb$omega_ka, pooled_tb$omega_ka_se)["pooled_est"]
        estimates[3, "rses"] <- rubin_fn(pooled_tb$omega_ka, pooled_tb$omega_ka_se)["pooled_rse"]
        estimates[4, "estimate"] <- rubin_fn(pooled_tb$omega_Cl, pooled_tb$omega_Cl_se)["pooled_est"]
        estimates[4, "rses"] <- rubin_fn(pooled_tb$omega_Cl, pooled_tb$omega_Cl_se)["pooled_rse"]
        estimates[5, "estimate"] <- rubin_fn(pooled_tb$b, pooled_tb$b_se)["pooled_est"]
        estimates[5, "rses"] <- rubin_fn(pooled_tb$b, pooled_tb$b_se)["pooled_rse"]
        
        # Save the results
        for (parami in params_for_estimation) {
          # Determine the action based on the suffix of `parami`
          action <- case_when(
            str_detect(parami, "_theta$") ~ "theta",
            str_detect(parami, "_theta_rse") ~ "theta_rse",
            str_detect(parami, "_eta$") ~ "eta",
            str_detect(parami, "_eta_rse") ~ "eta_rse",
            str_detect(parami, "b$") ~ "b",
            str_detect(parami, "b_rse") ~ "b_rse",
            TRUE ~ NA_character_ # Default case
          )
          
          if (!is.na(action)) {
            results_matrix[j, parami] <- switch(
              action,
              "theta" = pull(filter(estimates, param == gsub("_theta$", "", parami)), estimate),
              "theta_rse" = pull(filter(estimates, param == gsub("_theta_rse", "", parami)), rses),
              "eta" = pull(filter(estimates, param == paste0("omega_", gsub("_eta$", "", parami))), estimate),
              "eta_rse" = pull(filter(estimates, param == paste0("omega_", gsub("_eta_rse", "", parami))), rses),
              "b" = pull(filter(estimates, param == parami), estimate),
              "b_rse" = pull(filter(estimates, param == gsub("_rse", "", parami)), rses)
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
      path_to_save <- paste0("missingness/", imputation_method, "_", missing_percentage, "_", mechanism)
      write.csv(results_matrix, paste0(path_to_save, ".csv"), row.names = FALSE, quote = FALSE)
      write.csv(time_tb, paste0(path_to_save, "_time.csv"), row.names = FALSE)
    }
  }
}

