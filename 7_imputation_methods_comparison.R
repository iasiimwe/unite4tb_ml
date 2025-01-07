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
# omega_V <- omega_squared_fn(18.2)
# gamma_V <- omega_squared_fn(11.7)

sim_date <-"23_12_2024"
n_datasets <- 100
effect_size <- 0.5

imputation_methods <- c("CCA", "Mode")
mechanisms <- c("MCAR", "MAR", "MNAR")
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
        # Get DV data
        dv_dat <- read_csv(paste0("sim_DV/sim_data", j, "_", sim_date, ".csv"), show_col_types = FALSE) %>%
          select(!SNP1:SNP9) # Remove complete SNP data
        
        # Get missing SNP data
        dat_path <- paste0("datasets/missing/", missing_percentage, "/dat_", mechanism, "_", j)

        # Convert to ped format
        system(paste0("C:/Users/iasiimww/Downloads/plink/plink --bfile ", dat_path, " --recode tab --out temp"))
        
        # Process ped and map_files
        ped <- tibble(fread("temp.ped", header = FALSE))
        map <- tibble(fread("temp.map", header = FALSE))
        snps <- pull(map, V2)
        ped_start <- ped %>% select(V1:V6) # This selects FID, IID, father's ID, mother's ID, sex and phenotype data
        snp_columns <- colnames(ped)[!colnames(ped) %in% colnames(ped_start)]
        ped_snps <- ped %>%
          select(all_of(snp_columns)) 
        
        # True covariates
        true_snps <- fread(paste0("datasets/true_covar/true_snps_", j, ".txt"), header = FALSE, sep = " ")$V1
        # Add true to the true SNPs
        snps[snps %in% true_snps] <- paste0("true_", snps[snps %in% true_snps])
        # Change the ped_snps colnames
        colnames(ped_snps) <- snps 
        
        # Retain only the true SNPs, change to 'additive' coding, arrange them, and rename to SNP1-9
        ped_snps <- ped_snps %>%
          select(contains("true")) %>%
          mutate_all(dominant_coding_fn)
        true_snps <- gsub("[^a-zA-Z0-9 ]", "", true_snps) # Remove special characters
        colnames(ped_snps) <- gsub("[^a-zA-Z0-9 ]|true", "", colnames(ped_snps))
        
        ped_snps <- ped_snps %>%
          select(any_of(true_snps)) # This allows the order to be the same as in 'true_snps'
        colnames(ped_snps) <- paste0("SNP", 1:9)
        
        # Add back ID
        ped_snps$ID2 <- ped_start$V1 
        
        # Get time of code execution start
        start <- Sys.time()
        
        # Impute the NAs here
        #---------------------
        if (imputation_method == "CCA") {
          ped_snps <- ped_snps[complete.cases(ped_snps), ]
          # Note that the number of complete cases is different from the missingness we calculated e.g. example below:
            # Missing at least one SNP: round(sum(!complete.cases(ped_snps))/nrow(ped_snps) * 100, 2) # 27%
            # SNP percent missingness: i.e. round(sum(is.na(ped_snps))/prod(dim(ped_snps)) * 100, 2) # 3.56%
          if(nrow(ped_snps) < 10) next # At least ten participants per analysis - this is because with complete case analysis, numbers go down very quickly
          if(0 %in% apply(ped_snps[, sapply(ped_snps, is.numeric)], 2, max)) next # Means a SNP has only wild-type and MONOLIX will raise an error
          } else if (imputation_method == "Mode") { # Impute with wild-type (0)
          ped_snps <- ped_snps %>%
            mutate_at(vars(contains("SNP")), mode_x_fn)
            if(0 %in% apply(ped_snps[, sapply(ped_snps, is.numeric)], 2, max)) next # Means a SNP has only wild-type and MONOLIX will raise an error
        } else {
          stop("No method provided")
        }
        #---------------------
        
        # Join to the dv data
        dat <- ped_snps %>%
          left_join(dv_dat) %>%
          relocate(SNP1:ID2, .after = WT) %>%
          select(-NTIME, -ID2) 
        write.csv(dat, "temp.csv", row.names = FALSE)
        
        initializeLixoftConnectors(software = "monolix", path = "/pub59/iasiimwe/Lixoft/MonolixSuite2024R1/", force = TRUE)
        newProject(data = list(dataFile = 'temp.csv',
                               headerTypes = c("id", "time", "evid", "occ", "amount", "observation", 
                                               "catcov", "regressor", rep("catcov", 9))),
                   modelFile = 'tb_base_vinnard.txt')
        
        # Change error model, distribution, variability and correlations
        setErrorModel(DV = "combined1")
        setIndividualParameterDistribution(ka = "logNormal", TVV = "logNormal", TVCl = "logNormal")
        setIndividualParameterVariability(id = list (ka = TRUE, TVV = FALSE, TVCl = TRUE),
                                          OCC = list (ka = FALSE, TVV = FALSE, TVCl = FALSE))
        # getIndividualParameterModel() to confirm changes
        
        # Add sex and SNPs as covariates
        setCovariateModel(TVCl = c(SNP1 = TRUE, SNP2 = TRUE, SNP3 = TRUE, SNP4 = TRUE, 
                                   SNP5 = TRUE, SNP6 = TRUE, SNP7 = TRUE, SNP8 = TRUE, SNP9 = TRUE))
        
        # getPopulationParameterInformation()
        setPopulationParameterInformation(ka_pop = list(initialValue = 1.31, method = "MLE"),
                                          TVV_pop = list(initialValue = 28.57, method = "FIXED"),
                                          TVCl_pop = list(initialValue = 3.52, method = "MLE"),
                                          # beta_TVCl_SEX_2 = list(initialValue = -0.4, method = "FIXED"),
                                          beta_TVCl_SNP1_1 = list(initialValue = effect_size, method = "FIXED"),
                                          beta_TVCl_SNP2_1 = list(initialValue = effect_size, method = "FIXED"),
                                          beta_TVCl_SNP3_1 = list(initialValue = effect_size, method = "FIXED"),
                                          beta_TVCl_SNP4_1 = list(initialValue = effect_size, method = "FIXED"),
                                          beta_TVCl_SNP5_1 = list(initialValue = effect_size, method = "FIXED"),
                                          beta_TVCl_SNP6_1 = list(initialValue = effect_size, method = "FIXED"),
                                          beta_TVCl_SNP7_1 = list(initialValue = effect_size, method = "FIXED"),
                                          beta_TVCl_SNP8_1 = list(initialValue = effect_size, method = "FIXED"),
                                          beta_TVCl_SNP9_1 = list(initialValue = effect_size, method = "FIXED"),
                                          # gamma_TVV = list(initialValue = gamma_V, method = "MLE"),
                                          omega_TVCl = list(initialValue = omega_Cl, method = "MLE"),
                                          # omega_TVV = list(initialValue = omega_V, method = "MLE"),
                                          omega_ka = list(initialValue = omega_ka, method = "MLE"),
                                          a = list(initialValue = 2.41, method = "FIXED"),
                                          b = list(initialValue = 0.22, method = "MLE"),
                                          c = list(initialValue = 1, method = "FIXED"))
        
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




# Performance metrics
imputation_methods <- c("CCA", "Mode")
mechanisms <- c("MCAR", "MAR", "MNAR")
missing_percentages <- c(5, 10, 20, 50)
params_for_estimation <- c("ka_theta", "ka_theta_rse", "ka_eta", "ka_eta_rse", 
                           "Cl_theta", "Cl_theta_rse", "Cl_eta", "Cl_eta_rse", "b", "b_rse")

# Performance function
performance_fn <- function(parami, original_data = original_dat, imputed_data = imputed_dat) {
  # Successful run; examples of failures are when some SNPs become redundant 
  Success <- sum(rowSums(imputed_data == 0) != ncol(imputed_data)) * 100/nrow(imputed_data)
  
  # Replace rows with all zeros in the imputed data to NA
  imputed_data[rowSums(imputed_data == 0) == ncol(imputed_data), ] <- NA
  
  # Get relevant columns
  true <- original_data[, parami]
  estimated <- imputed_data[, parami]
  
  # Remove the NAs
  true <- true[!is.na(estimated)]
  estimated <- estimated[!is.na(estimated)]
  
  if (sum(is.na(estimated)) != length(estimated)) {
    # Bias (For acceptable performance we use an upper limit for PB of 5%. (Demirtas, Freels, and Yucel 2008))
    eMLAR <- (exp(mean(log(estimated / true))) - 1) * 100 # Handles proportional biases well, especially for log-normally distributed data.
    MRPE <- mean((estimated - true) / true) * 100 # Evaluates bias for individual observations.
    rMPE <- (mean(estimated - true) / mean(true)) * 100 # Evaluates overall bias in mean predictions vs. mean true.
    
    # Precision
    eMALAR <- (exp(mean(abs(log(estimated / true)))) - 1) * 100 # Evaluates proportional precision, especially for log-normally distributed data.
    RMSRE <- (sqrt(mean(((estimated - true) / true)^2))) * 100 # Penalizes large deviations and outliers more heavily.
    MAPE <- mean(abs((estimated - true) / true)) * 100 # Measures average magnitude of relative errors (robust to outliers).

      } else { # All rows are NA
    eMLAR <- NA; MRPE <- NA; rMPE <- NA; eMALAR <- NA; RMSRE <- NA; MAPE <- NA
  }
  return(tibble(data.frame(cbind(Success, eMLAR, MRPE, rMPE, eMALAR, RMSRE, MAPE))))
}

# Additional metrics (https://stefvanbuuren.name/fimd/sec-evaluation.html)
  # Refers to Raw Bias (RB) as mean(estimated - true), and Percent Bias (PB) as MAPE
  # Coverage rate (CR) - the proportion of confidence intervals that contain the true value.
     # CR <- rowMeans(res[,, "2.5 %"] < true & true < res[,, "97.5 %"])
  # Average width (AW) - the average width of the confidence interval, an indicator of statistical efficiency.
     # AW <- rowMeans(res[,, "97.5 %"] - res[,, "2.5 %"])
performance_fn_mice <- function(parami, original_data = original_dat, imputed_data = imputed_dat) {
  # Replace rows with all zeros in the imputed data to NA
  imputed_data[rowSums(imputed_data == 0) == ncol(imputed_data), ] <- NA
  
  # Get relevant columns
  true <- pull(original_data[, parami])
  estimated <- imputed_data[, parami]
  estimated_se <- original_data[, paste0(parami, "_rse")]
  estimated <- bind_cols(estimated, estimated_se)
  colnames(estimated) <- c("paramj", "paramj_rse")
  estimated <- estimated %>%
    mutate(se = paramj * paramj_rse/ 100,
           lower_ci = paramj + (qnorm(0.025) * se),
           upper_ci = paramj + (qnorm(0.975) * se),
           truej = true,
           AW = upper_ci - lower_ci,
           CR = if_else(lower_ci <= truej & upper_ci >= truej, 1, 0)) %>%
    filter(!is.na(paramj)) # Remove the NAs
  
  if (nrow(estimated) > 0) {
    AW <- mean(estimated$AW); CR <- mean(estimated$CR) # Don't multiply CR by 100, so that we can plot them together with AW
    } else {
      AW <- NA; CR <- NA
    }
  return(tibble(data.frame(cbind(AW, CR))))
} 

# Get original data results
original_dat <- read_csv(paste0("missingness/complete.csv"), show_col_types = FALSE)

results <- tibble()
for (imputation_method in imputation_methods) {
  for (mechanism in mechanisms) {
    for (missing_percentage in missing_percentages) {
      
      # Get imputed data results
      path_to_save <- paste0("missingness/", imputation_method, "_", missing_percentage, "_", mechanism)
      imputed_dat <- read_csv(paste0(path_to_save, ".csv"), show_col_types = FALSE)
      
      for (param_for_estimation in params_for_estimation) {
        
        resultsi <- performance_fn(param_for_estimation) %>%
          mutate(Analysis = paste0(c(imputation_method, mechanism, missing_percentage, param_for_estimation), collapse = "|"))
        if (!str_detect(param_for_estimation, "rse")) {
          resultsj <- performance_fn_mice(param_for_estimation)
          resultsi <- bind_cols(resultsi, resultsj)
        }
        results <- bind_rows(results, resultsi)
      }
    }
  }
}

# Plots by missing mechanisms
for (mechanism in mechanisms) {
  resultsi <- results %>%
    filter(str_detect(Analysis, mechanism)) %>%
    mutate(Analysis = gsub(mechanism, "", Analysis),
           Analysis = gsub("\\|\\|", "|", Analysis)) %>%
    separate(Analysis, into = c("Method", "Percentage", "Parameter"), sep = "\\|") %>%
    mutate(Label = paste(Method, " (", Percentage, "%)", sep = ""),
           Parameter = gsub("_theta", " THETA", Parameter),
           Parameter = gsub("_eta", " ETA", Parameter),
           Parameter = gsub("_rse", " RSE", Parameter),
           Parameter = factor(Parameter, levels = c("Cl THETA", "Cl THETA RSE", "Cl ETA", "Cl ETA RSE",
                                                    "ka THETA", "ka THETA RSE", "ka ETA", "ka ETA RSE",
                                                    "b", "b RSE")),
           Label = factor(Label, levels = unique(Label)))
  
  # Reshape the data to a long format for plotting
  resultsi_long <- resultsi %>%
    pivot_longer(cols = c(Success, eMLAR, MRPE, rMPE, eMALAR, RMSRE, MAPE),
                 names_to = "Metric",
                 values_to = "Value") 
  
  # Create the barplots (start with success)
  resultsi_long %>%
    filter(Metric == "Success") %>%
    select(-Parameter) %>%
    distinct() %>%
    ggplot(aes(x = Label, y = Value)) +
    geom_bar(stat = "identity", position = "dodge", fill = "lightblue") +
    theme_bw() +
    labs(
      title = paste0("Barplots for Successful Runs (", mechanism, ")"),
      x = "Missing Data Mechanism (Missingness Percentage)",
      y = "Percentage"
    ) +
    theme(axis.text = element_text(size = 10, colour = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title = element_text(size = 12, face = "bold"),
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

  x <- 0.55
  ggsave(paste0("Success_", mechanism, ".png"), height = 9 * x, width = 16 * x)
  
  # Bias metrics
  resultsi_long %>%
    filter(Metric %in% c("eMLAR", "MRPE", "rMPE")) %>%
    ggplot(aes(x = Label, y = Value, fill = Metric)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~ Parameter, scales = "free_y", nrow = 3) +
    theme_bw() +
    labs(
      title = paste0("Barplots of Bias Metrics by Parameter and Method (", mechanism, ")"),
      x = "Missing Data Mechanism (Missingness Percentage)",
      y = "Percentage"
    ) +
    theme(legend.position = c(0.9, 0.1),            # Set legend position to bottom-right
        #legend.justification = c(1, 0),          # Adjust the anchor point of the legend
        legend.background = element_rect(fill = "white", colour = "black"), # Optional: Add a white background to the legend
        strip.text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10, colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  x <- 1
  ggsave(paste0("Bias_", mechanism, ".png"), height = 9 * x, width = 16 * x)
  
  # Precision metrics
  resultsi_long %>%
    filter(Metric %in% c("eMALAR", "MAPE", "RMSRE")) %>%
    ggplot(aes(x = Label, y = Value, fill = Metric)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~ Parameter, scales = "free_y", nrow = 3) +
    theme_bw() +
    labs(
      title = paste0("Barplots of Precision Metrics by Parameter and Method (", mechanism, ")"),
      x = "Missing Data Mechanism (Missingness Percentage)",
      y = "Percentage"
    ) +
    theme(legend.position = c(0.9, 0.1),            # Set legend position to bottom-right
          #legend.justification = c(1, 0),          # Adjust the anchor point of the legend
          legend.background = element_rect(fill = "white", colour = "black"), # Optional: Add a white background to the legend
          strip.text = element_text(size = 12, face = "bold"),
          axis.text = element_text(size = 10, colour = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title = element_text(size = 12, face = "bold"),
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  x <- 1
  ggsave(paste0("Precision_", mechanism, ".png"), height = 9 * x, width = 16 * x)
}


# Combined plots (no RSEs)
# ------------------------
resultsi_long <- results %>%
  mutate(Analysis = gsub("\\|\\|", "|", Analysis)) %>%
  separate(Analysis, into = c("Method", "Mechanism", "Percentage", "Parameter"), sep = "\\|") %>%
  mutate(Label = paste(Method, " (", Percentage, "%)", sep = ""),
         Parameter = gsub("_theta", " THETA", Parameter),
         Parameter = gsub("_eta", " ETA", Parameter),
         Parameter = gsub("_rse", " RSE", Parameter),
         Parameter = factor(Parameter, levels = c("Cl THETA", "Cl THETA RSE", "Cl ETA", "Cl ETA RSE",
                                                  "ka THETA", "ka THETA RSE", "ka ETA", "ka ETA RSE",
                                                  "b", "b RSE")),
         Label = factor(Label, levels = unique(Label)),
         Mechanism = factor(Mechanism, levels = unique(Mechanism))) %>%
  pivot_longer(cols = c(Success, eMLAR, MRPE, rMPE, eMALAR, RMSRE, MAPE, AW, CR),
               names_to = "Metric",
               values_to = "Value")

# Success
resultsi_long %>%
  filter(Metric == "Success") %>%
  select(-Parameter) %>%
  distinct() %>%
  ggplot(aes(x = Label, y = Value)) +
  geom_bar(stat = "identity", position = "dodge", fill = "lightblue") +
  # Add labels below bars for values < 100 and > 0
  geom_text(
    aes(label = ifelse(Value < 100 & Value > 0, round(Value, 1), "")),
    position = position_dodge(width = 0.9),
    vjust = 1.5,  # Position the label below the bar
    size = 3.5
  ) +
  facet_wrap(~ Mechanism, scales = "free_y", nrow = 2) +
  theme_bw() +
  labs(
    title = "Barplots for Successful Runs",
    x = "Imputation Method (Missingness Percentage)",
    y = "Percentage"
  ) +
  theme(strip.text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10, colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
x <- 0.65
ggsave("Success.png", height = 9 * x, width = 16 * x)

# Bias (no RSEs)
resultsi_long %>%
  filter(Metric %in% c("eMLAR", "MRPE", "rMPE")) %>%
  filter(!str_detect(Parameter, "RSE")) %>%
  mutate(Parameter = factor(Parameter, levels = c("Cl THETA", "Cl ETA",
                                                  "ka THETA", "ka ETA", "b"))) %>%
  ggplot(aes(x = Label, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_nested_wrap(vars(Mechanism, Parameter), nrow = 3, scales = "free_y") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  theme_bw() +
  labs(
    # title = "Barplots of Bias Metrics by Parameter and Method",
    x = "Imputation Method (Missingness Percentage)",
    y = "Percentage"
  ) +
  theme(legend.position = "top",            # Set legend position to bottom-right
        #legend.justification = c(1, 0),          # Adjust the anchor point of the legend
        legend.background = element_rect(fill = "white", colour = "black"), # Optional: Add a white background to the legend
        strip.text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10, colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 12, face = "bold")
        # ,
        # plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
        )
x <- 1
ggsave("Bias.png", height = 9 * x, width = 16 * x)
  
# Precision (no RSEs)
resultsi_long %>%
  filter(Metric %in% c("eMALAR", "MAPE", "RMSRE")) %>%
  filter(!str_detect(Parameter, "RSE")) %>%
  mutate(Parameter = factor(Parameter, levels = c("Cl THETA", "Cl ETA",
                                                  "ka THETA", "ka ETA", "b"))) %>%
  ggplot(aes(x = Label, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_nested_wrap(vars(Mechanism, Parameter), nrow = 3, scales = "free_y") +
  # geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  theme_bw() +
  labs(x = "Imputation Method (Missingness Percentage)", y = "Percentage") +
  theme(legend.position = "top",           
        legend.background = element_rect(fill = "white", colour = "black"), # Optional: Add a white background to the legend
        strip.text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10, colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 12, face = "bold"))
x <- 1
ggsave("Precision.png", height = 9 * x, width = 16 * x)


# Additional metrics (AW and CR)
resultsi_long %>%
  filter(Metric %in% c("AW", "CR")) %>%
  filter(!str_detect(Parameter, "RSE")) %>%
  mutate(Parameter = factor(Parameter, levels = c("Cl THETA", "Cl ETA",
                                                  "ka THETA", "ka ETA", "b")),
         Metric = case_when(Metric == "AW" ~ "Average width", Metric == "CR" ~ "Coverage rate" )) %>%
  ggplot(aes(x = Label, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_nested_wrap(vars(Mechanism, Parameter), nrow = 3, scales = "free_y") +
  geom_hline(yintercept = 0.90, linetype = "dashed", color = "black") +
  theme_bw() +
  labs(x = "Imputation Method (Missingness Percentage)", y = "Percentage") +
  theme(legend.position = "top",           
        legend.background = element_rect(fill = "white", colour = "black"), # Optional: Add a white background to the legend
        strip.text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10, colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 12, face = "bold"))
x <- 1
ggsave("AW_CR.png", height = 9 * x, width = 16 * x)
