# Load required packages
# ---------------------
library(lixoftConnectors)
initializeLixoftConnectors(software = "monolix", path = "/pub59/iasiimwe/Lixoft/MonolixSuite2024R1/", force = TRUE)
library(tidyverse)
library(data.table)
library(ggh4x)

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

# Imputing using the mode/wildtype
mode_x_fn <- function(x) {
  x[is.na(x)] <- 0
  return(x)
}

## Random effects (considering log distribution)
omega_squared_fn <- function(percent_CV) log((percent_CV/100)^2 + 1)

# # Without re-inflation
omega_ka <- omega_squared_fn(68.5)
omega_Cl <- omega_squared_fn(17.5)

sim_date <-"23_12_2024"
n_datasets <- 100
effect_size <- 0.5

imputation_methods <- c("CCA", "Mode", "michigan_afr", "michigan_eur")
mechanisms <- c("MCAR", "MAR", "MNAR")
missing_percentages <- c(5, 10, 20, 50)

for (k in seq_along(imputation_methods)) {
  imputation_method <- imputation_methods[[k]]
  
  for (mechanism in mechanisms) {
    for (missing_percentage in missing_percentages) {
      if (all(imputation_method %in% c("michigan_afr", "michigan_eur") & mechanism %in% c("MAR", "MNAR"))) next
      if (all(imputation_method %in% c("michigan_afr", "michigan_eur") & missing_percentages %in% c(10, 20, 50))) next
      
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
        # Get DV data
        dv_dat <- read_csv(paste0("sim_DV/sim_data", j, "_", sim_date, ".csv"), show_col_types = FALSE) %>%
          select(!SNP1:SNP9) # Remove complete SNP data
        
        # Get missing SNP data
        dat_path <- paste0("missing/", missing_percentage, "/dat_", mechanism, "_", j)

        # Convert to ped format
        system(paste0("/pub59/iasiimwe/plink1.9/plink --bfile ", dat_path, " --recode tab --out temp_comparison"))
        
        # Process ped and map_files
        ped <- tibble(fread("temp_comparison.ped", header = FALSE))
        map <- tibble(fread("temp_comparison.map", header = FALSE))
        snps <- pull(map, V2)
        ped_start <- ped %>% select(V1:V6) # This selects FID, IID, father's ID, mother's ID, sex and phenotype data
        snp_columns <- colnames(ped)[!colnames(ped) %in% colnames(ped_start)]
        ped_snps <- ped %>%
          select(all_of(snp_columns)) 
        colnames(ped_snps) <- snps 
        
        # True covariates
        true_snps <- fread(paste0("true_covar/true_snps_", j, ".txt"), header = FALSE, sep = " ")$V1

        # Retain only the true SNPs, change to 'dominant' coding, arrange them
        ped_snps <- ped_snps %>%
          select(any_of(true_snps)) %>% # This allows the order to be the same as in 'true_snps'
          mutate_all(dominant_coding_fn) %>%
          mutate(ID2 = ped_start$V1)

        # Get time of code execution start
        start <- Sys.time()
        
        # Impute the NAs here
        #---------------------
        if (imputation_method == "CCA") {
          # Rename the SNPs to SNP1-9
          colnames(ped_snps) <- c(paste0("SNP", 1:9), "ID2")
          ped_snps <- ped_snps[complete.cases(ped_snps), ]
          # Note that the number of complete cases is different from the missingness we calculated e.g. example below:
          # Missing at least one SNP: round(sum(!complete.cases(ped_snps))/nrow(ped_snps) * 100, 2) # 27%
          # SNP percent missingness: i.e. round(sum(is.na(ped_snps))/prod(dim(ped_snps)) * 100, 2) # 3.56%
          if(nrow(ped_snps) < 10) next # At least ten participants per analysis - this is because with complete case analysis, numbers go down very quickly
          # if(0 %in% apply(ped_snps[, sapply(ped_snps, is.numeric)], 2, max)) next # Means a SNP has only wild-type and MONOLIX will raise an error
        } else if (imputation_method == "Mode") { # Impute with wild-type (0)
          # Rename the SNPs to SNP1-9
          colnames(ped_snps) <- c(paste0("SNP", 1:9), "ID2")
          ped_snps <- ped_snps %>%
            mutate_at(vars(contains("SNP")), mode_x_fn)
          # if(0 %in% apply(ped_snps[, sapply(ped_snps, is.numeric)], 2, max)) next # Means a SNP has only wild-type and MONOLIX will raise an error
        } else if (imputation_method %in% c("michigan_afr", "michigan_eur")) {
          ref_population <- gsub("michigan_", "", imputation_method)
          # Let us get ids
          job_id <- read_csv(paste0("michigan/job_id_", ref_population, ".csv"), show_col_types = FALSE)$job_id[[j]]
          # Check if the id was successful (one reason for failure)
          job_id_path <- paste0("/pub59/iasiimwe/imputationserver/workspace/", job_id, "/output/chr_1.zip")
          if(!file.exists(job_id_path)) next
          
          # Get true SNPs
          true_snps_m <- paste0("SNP_", fread(paste0("true_covar/true_snps_", j, "_GRCh37.txt"), header = TRUE, sep = " ")$GRCh37_BP)
          
          # Process ped and map_files
          ped_path <- paste0("michigan/chr", j, "_", job_id, "_extracted")
          ped_m <- tibble(fread(paste0(ped_path, ".ped"), header = FALSE))
          map_m <- tibble(fread(paste0(ped_path, ".map"), header = FALSE))
          snps_m <- paste0("SNP_", pull(map_m, V4))
          ped_start_m <- ped_m %>% select(V1:V6) # This selects FID, IID, father's ID, mother's ID, sex and phenotype data
          snp_columns_m <- colnames(ped_m)[!colnames(ped_m) %in% colnames(ped_start_m)]
          
          ped_snps_m <- ped_m %>%
            select(all_of(snp_columns_m)) 
          colnames(ped_snps_m) <- snps_m 
          
          # Select the true SNPs while re-ordering them
          ped_snps_m <- ped_snps_m %>%
            select(any_of(true_snps_m)) 
          
          # Check if all 9 SNPs were imputed (another reason for failure - track both)
          if(ncol(ped_snps_m) < 9) next
          
          # Dominant codings
          ped_snps <- ped_snps_m %>%
            mutate_all(dominant_coding_fn)  %>%
            mutate(ID2 = ped_start_m$V1)
          # if(0 %in% apply(ped_snps[, sapply(ped_snps, is.numeric)], 2, max)) next # Means a SNP has only wild-type and MONOLIX will raise an error
          ped_snps_m <- NULL
          colnames(ped_snps) <- c(paste0("SNP", 1:9), "ID2")
        } else {
          stop("MICE methods")
        }
        #---------------------
        # Find columns that have only wild-type
        snps_to_remove <- colnames(ped_snps)[apply(ped_snps[, sapply(ped_snps, is.numeric)], 2, max) == 0]
        snps_to_remove <- snps_to_remove[str_detect(snps_to_remove, "SNP")] # Not to remove ID
        
        # Join to the dv data
        dat <- ped_snps %>%
          select(!all_of(snps_to_remove)) %>%
          left_join(dv_dat) %>%
          relocate(ID:WT, .before = 1) %>%
          select(-NTIME, -ID2) 
        write.csv(dat, "temp_monolix.csv", row.names = FALSE)
        
        initializeLixoftConnectors(software = "monolix", path = "/pub59/iasiimwe/Lixoft/MonolixSuite2024R1/", force = TRUE)
        newProject(data = list(dataFile = 'temp_monolix.csv',
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
        estimates <- tibble(param = getEstimatedStandardErrors()$stochasticApproximation$parameter,
                            estimate = getEstimatedPopulationParameters()[names(getEstimatedPopulationParameters()) %in% getEstimatedStandardErrors()$stochasticApproximation$parameter],
                            rses = getEstimatedStandardErrors()$stochasticApproximation$rse) %>%
          mutate(param = gsub("TV|_pop", "", param))
        
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

