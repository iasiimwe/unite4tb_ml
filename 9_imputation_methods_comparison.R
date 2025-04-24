# Load required packages
# ---------------------
library(tidyverse)
library(data.table)
library(ggh4x)
library(ggbreak)
library(RColorBrewer)
library(patchwork)
library(export)

select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate

# Performance function
performance_fn <- function(parami, original_data = original_dat, imputed_data = imputed_dat) {
  # Successful run; examples of failures are when some SNPs become redundant 
  imputed_data <- imputed_data %>% select(!contains("df")) # To account for the Infs added to the single imputation methods
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
  # Degrees of freedom
  df_data <- imputed_data %>% select(contains("df"))
  imputed_data <- imputed_data %>% select(!contains("df")) 
  df_data[rowSums(imputed_data == 0) == ncol(imputed_data), ] <- NA
  
  # Replace rows with all zeros in the imputed data to NA
  imputed_data[rowSums(imputed_data == 0) == ncol(imputed_data), ] <- NA
  
  # Get relevant columns
  true <- pull(original_data[, parami])
  estimated <- imputed_data %>% bind_cols(df_data) %>% select(contains(parami))
  colnames(estimated) <- c("paramj", "paramj_rse", "paramj_df")
  estimated <- estimated %>%
    mutate(se = paramj * paramj_rse/ 100,
           lower_ci = paramj + (qt(0.025, paramj_df) * se), 
           upper_ci = paramj + (qt(0.975, paramj_df) * se),
           truej = true,
           AW = upper_ci - lower_ci,
           CR = if_else(lower_ci <= truej & upper_ci >= truej, 1, 0)) %>%
    filter(!is.na(paramj)) # Remove the NAs

  if (nrow(estimated) > 0) {
    AW <- mean(estimated$AW); CR <- mean(estimated$CR) * 100
  } else {
    AW <- NA; CR <- NA
  }
  return(tibble(data.frame(cbind(AW, CR))))
} 

# Time function
time_fn <- function(imputation_method = NULL, mechanism = NULL, missing_percentage = NULL, 
                    original = TRUE) {
  # Working environment assumed to have relevant the folders
  if (original) {
    time_tb <- read_csv(paste0("missingness_", effect_scenario, "/original_data_time.csv"), show_col_types = FALSE) %>%
      mutate(end = ymd_hms(gsub("T ", "", end)),
             start = ymd_hms(gsub("T ", "", start)),
             time_s = as.numeric(end-start, "secs")) %>%
      select(time_s)
    return(time_tb)
  } else {
    # All methods (Monolix fitting)
    path_to_imputed_data <- ifelse(imputation_method %in% c("michigan_afr_single"),
                                   paste0("missingness_", effect_scenario, "/", gsub("_single", "", imputation_method), "_", 
                                          missing_percentage, "_", mechanism, "_single_time.csv"),
                                   paste0("missingness_", effect_scenario, "/", gsub("_0.*", "", imputation_method), "_", missing_percentage, "_", mechanism, "_time.csv"))
    time_tb <- read_csv(path_to_imputed_data, show_col_types = FALSE) %>%
      mutate(end = ymd_hms(gsub("T ", "", end)),
             start = ymd_hms(gsub("T ", "", start)),
             time_s = as.numeric(end-start, "secs")) %>%
      select(time_s)
    
    # Imputation server (single, multiple with pre-filtering)
    if (imputation_method %in% c("michigan_afr_single", "michigan_afr_ms")) {
      # Minimac4 preparation
      time_tb <- read_csv("michigan/time.csv", show_col_types = FALSE) %>%
        mutate(end = ymd_hms(gsub("T ", "", end)),
               start = ymd_hms(gsub("T ", "", start)),
               time = as.numeric(end-start, "secs")) %>%
        select(time) %>%
        bind_cols(time_tb) %>%
        mutate(time_s = time_s + time) %>%
        select(time_s)

      # Minimac4 imputation server
      ref_population <- gsub("_single|_ms", "", gsub("michigan_", "", imputation_method))
      time_tb2 <- rep(NA_real_, n_datasets)
      for (i in seq_along(time_tb2)) {
        job_id <- read_csv(paste0("michigan/job_id_", ref_population, ".csv"), show_col_types = FALSE)$job_id[[i]]
        job_id_path <- paste0("/pub59/iasiimwe/imputationserver/workspace/", job_id, "/output/chr_1.zip")
        if(!file.exists(job_id_path)) next
        job_txt <- readLines(paste0("/pub59/iasiimwe/imputationserver/workspace/", job_id, "/logs/job.txt"))
        time_tb2[[i]] <- as.numeric(ymd_hms(gsub(" Cleanup.*", "", job_txt[length(job_txt)])) - ymd_hms(gsub(" Setup.*", "", job_txt[1])), "secs")
      }
      time_tb <- time_tb %>%
        mutate(time_s = time_s + time_tb2)

      # Minimac4 results extraction
      time_tb <- tibble(dataset = 1:n_datasets) %>%
        left_join(read_csv(paste0("michigan/time_", ref_population, ".csv"), show_col_types = FALSE)) %>%
        mutate(time = as.numeric(ymd_hms(end) - ymd_hms(start), "secs")) %>%
        select(time) %>%
        bind_cols(time_tb) %>%
        mutate(time_s = time_s + time) %>%
        select(time_s)
      
      if (imputation_method == "michigan_afr_ms") {
        # Info score extraction
        time_tb <- tibble(dataset = 1:n_datasets) %>%
          left_join(read_csv(paste0("michigan/info_scores_time_", ref_population, ".csv"), show_col_types = FALSE)) %>%
          mutate(time = as.numeric(ymd_hms(end) - ymd_hms(start), "secs")) %>%
          select(time) %>%
          bind_cols(time_tb) %>%
          mutate(time_s = time_s + time) %>%
          select(time_s) 
      }
    }
    
    # Imputation server (multiple)
    if (imputation_method %in% c("michigan_afr", "michigan_afr_5dat", "michigan_afr_5dat")) {
        ref_population <- gsub("_\\d.*", "", gsub("michigan_", "", imputation_method))
      # Minimac4 preparation
      time_tb <- read_csv(paste0("michigan/", mechanism, "_", missing_percentage, "/time.csv"), show_col_types = FALSE)  %>%
        mutate(end = ymd_hms(gsub("T ", "", end)),
               start = ymd_hms(gsub("T ", "", start)),
               time = as.numeric(end-start, "secs")) %>%
        select(time) %>%
        bind_cols(time_tb) %>%
        mutate(time_s = time_s + time) %>%
        select(time_s)

      # Minimac4 imputation server
      time_tb2 <- rep(NA_real_, n_datasets)
      for (i in seq_along(time_tb2)) {
        job_id <- read_csv(paste0("michigan/job_id_", mechanism, "_", missing_percentage, "_", ref_population, ".csv"), show_col_types = FALSE)$job_id[[i]]
        job_id_path <- paste0("/pub59/iasiimwe/imputationserver/workspace/", job_id, "/output/chr_1.zip")
        if(!file.exists(job_id_path)) next
        job_txt <- readLines(paste0("/pub59/iasiimwe/imputationserver/workspace/", job_id, "/logs/job.txt"))
        time_tb2[[i]] <- as.numeric(ymd_hms(gsub(" Cleanup.*", "", job_txt[length(job_txt)])) - ymd_hms(gsub(" Setup.*", "", job_txt[1])), "secs")
      }
      time_tb <- time_tb %>%
        mutate(time_s = time_s + time_tb2)
      
      # Minimac4 results extraction
      time_tb <- tibble(dataset = 1:n_datasets) %>%
        left_join(read_csv(paste0("michigan/time_", mechanism, "_", missing_percentage, "_", ref_population, ".csv"), show_col_types = FALSE)) %>%
        mutate(time = as.numeric(ymd_hms(end) - ymd_hms(start), "secs")) %>%
        select(time) %>%
        bind_cols(time_tb) %>%
        mutate(time_s = time_s + time) %>%
        select(time_s)
      
      # Info score extraction
      time_tb <- tibble(dataset = 1:n_datasets) %>%
        left_join(read_csv(paste0("michigan/info_scores_time_", ref_population, "_", mechanism, "_", missing_percentage, ".csv"), show_col_types = FALSE)) %>%
        mutate(time = as.numeric(ymd_hms(end) - ymd_hms(start), "secs")) %>%
        select(time) %>%
        bind_cols(time_tb) %>%
        mutate(time_s = time_s + time) %>%
        select(time_s)
    }
    
    # Mice imputation
    if (imputation_method %in% c("pmm", "midastouch", "sample", "cart", "rf", "mean", "norm", "ri")) {
      path_to_mice <- paste0("mice/", imputation_method, "_", missing_percentage, "_", mechanism)
      time_tb2 <- if(mechanism != "MNAR_post") {
        read_csv(paste0(path_to_mice, "_time.csv"), show_col_types = FALSE) 
        } else read_csv(paste0(gsub("_post", "", path_to_mice), "_time_post.csv"), show_col_types = FALSE)
      time_tb <- time_tb2 %>%
        mutate(end = ymd_hms(gsub("T ", "", end)),
               start = ymd_hms(gsub("T ", "", start)),
               time = as.numeric(end-start, "secs")) %>%
        select(time) %>%
        bind_cols(time_tb) %>%
        mutate(time_s = time_s + time) %>%
        select(time_s)
    }
    
    # Chat GPT
    if (imputation_method %in% c("chatgpt")) {
      path_to_chatgpt <- paste0("chatgpt/", imputation_method, "_", missing_percentage, "_", mechanism)
      time_tb2 <- read_csv(paste0(path_to_mice, "_time.csv"), show_col_types = FALSE) 
      time_tb <- time_tb2 %>%
        mutate(end = ymd_hms(gsub("T ", "", end)),
               start = ymd_hms(gsub("T ", "", start)),
               time = as.numeric(end-start, "secs")) %>%
        select(time) %>%
        bind_cols(time_tb) %>%
        mutate(time_s = time_s + time) %>%
        select(time_s)
    }
    
    return(time_tb)
  }
}

n_datasets <- 100
for (effect_scenario in c("low", "high")) {
  # Get original data results
  original_dat <- read_csv(paste0("missingness_", effect_scenario, "/complete.csv"), show_col_types = FALSE)
  original_time <- time_fn()
  
  # Performance metrics
  imputation_methods <- c("CCA", "Mode", "michigan_afr_single", "michigan_afr_ms", "michigan_afr", 
                          # "michigan_afr_5dat", "michigan_afr_50dat",
                          "pmm", "midastouch", "sample", "cart", "rf", "mean", "norm", "ri")
  mechanisms <- c("MCAR", "MAR", "MNAR")
  missing_percentages <- c(5, 10, 20, 50)
  params_for_estimation <- c("ka_theta", "ka_theta_rse", "ka_eta", "ka_eta_rse", 
                             "Cl_theta", "Cl_theta_rse", "Cl_eta", "Cl_eta_rse", "b", "b_rse")
  results <- tibble()
  results_time <- tibble()
  for (imputation_method in imputation_methods) {
    for (mechanism in mechanisms) {
      for (missing_percentage in missing_percentages) {
        if (all(imputation_method %in% c("michigan_afr_single", "michigan_afr_ms") & mechanism %in% c("MAR", "MNAR"))) next
        if (all(imputation_method %in% c("michigan_afr_single") & missing_percentage %in% c(10, 20, 50))) next
        if (all(imputation_method %in% c("michigan_afr", "michigan_afr_5dat", "michigan_afr_50dat", "michigan_afr_ms") & missing_percentage %in% c(20, 50))) next
        
        message(paste0("Starting\n     Method: ",imputation_method, "\n     Mechanism: ", 
                       mechanism, "\n     missing %: ", missing_percentage))
        
        # Get imputed data results
        path_to_imputed_data <- ifelse(imputation_method %in% c("michigan_afr_single"),
                                       paste0("missingness_", effect_scenario, "/", gsub("_single", "", imputation_method), "_", 
                                              missing_percentage, "_", mechanism, "_single"),
                                       paste0("missingness_", effect_scenario, "/", gsub("_0.*", "", imputation_method), "_", missing_percentage, "_", mechanism))
        imputed_dat <- read_csv(paste0(path_to_imputed_data, ".csv"), show_col_types = FALSE)
        
        if (imputation_method %in% c("CCA", "Mode", "michigan_afr_single")) {
          dfs = 100 - 5 # df = n - k, where k = 5 estimated parameters (CL_theta, CL_eta, ka_theta, ka_eta, b)
          # This is sligthly more accurate than using df = Inf (z-distribution for singly imputed data, n = 100)
          imputed_dat <- imputed_dat %>%
            mutate(ka_theta_df = dfs, .after = ka_theta_rse) %>% 
            mutate(ka_eta_df = dfs, .after = ka_eta_rse) %>%
            mutate(Cl_theta_df = dfs, .after = Cl_theta_rse) %>%
            mutate(Cl_eta_df = dfs, .after = Cl_eta_rse) %>%
            mutate(b_df = dfs, .after = b_rse)
        }
        
        for (param_for_estimation in params_for_estimation) {
          resultsi <- performance_fn(param_for_estimation) %>%
            mutate(Analysis = paste0(c(imputation_method, mechanism, missing_percentage, param_for_estimation), collapse = "|"))
          if (!str_detect(param_for_estimation, "rse")) {
            resultsj <- performance_fn_mice(param_for_estimation)
            resultsi <- bind_cols(resultsi, resultsj)
          }
          results <- bind_rows(results, resultsi)
        }
        
        # Get time results
        results_timei <- time_fn(imputation_method, mechanism, missing_percentage, original = FALSE) %>%
          mutate(Analysis = paste0(c(imputation_method, mechanism, missing_percentage), collapse = "|"))
        results_time <- bind_rows(results_time, results_timei)
        
        message("Completed!")
      }
    }
  }
  
  # Change to long format and save the results
  results <- results %>%
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
  results_time <- results_time %>%
    mutate(Analysis = gsub("\\|\\|", "|", Analysis)) %>%
    separate(Analysis, into = c("Method", "Mechanism", "Percentage"), sep = "\\|") %>%
    bind_rows(mutate(original_time, Method = "Reference", Mechanism = "MCAR", Percentage = "5")) %>% # Place holders
    mutate_if(is.character, factor)
  write_csv(results, paste0("performance_metrics_", effect_scenario, ".csv"))
  write_csv(results_time, paste0("performance_time_", effect_scenario, ".csv"))
}

# Get detailed results for multiple imputation michigan MCAR afr (10 datasets) high effect scenario
# --------------------------------------------------------------------------------------------------
effect_scenario <- "low"
imputation_method <- "michigan_afr"
mechanism <- "MCAR"
missing_percentage <- 10
n_datasets <- 100

# Monolix fitting
path_to_imputed_data <- paste0("missingness_", effect_scenario, "/", gsub("_0.*", "", imputation_method), "_", missing_percentage, "_", mechanism, "_time.csv")
monolix_fitting <- read_csv(path_to_imputed_data, show_col_types = FALSE) %>%
  mutate(end = ymd_hms(gsub("T ", "", end)),
         start = ymd_hms(gsub("T ", "", start)),
         time_s = as.numeric(end-start, "secs")) %>%
  select(time_s) %>%
  mutate(dataset = row_number(),
         Step = "Model fitting (Monolix)")

ref_population <- gsub("_0.*", "", gsub("michigan_", "", imputation_method))
# Minimac4 preparation
minimac <- read_csv(paste0("michigan/", mechanism, "_", missing_percentage, "/time.csv"), show_col_types = FALSE)  %>%
  mutate(end = ymd_hms(gsub("T ", "", end)),
         start = ymd_hms(gsub("T ", "", start)),
         time_s = as.numeric(end-start, "secs")) %>%
  select(time_s) %>%
  mutate(dataset = row_number(),
         Step = "Preparation for imputation")

# Minimac4 imputation server
time_tb2 <- rep(NA_real_, n_datasets)
for (i in seq_along(time_tb2)) {
  job_id <- read_csv(paste0("michigan/job_id_", mechanism, "_", missing_percentage, "_", ref_population, ".csv"), show_col_types = FALSE)$job_id[[i]]
  job_id_path <- paste0("/pub59/iasiimwe/imputationserver/workspace/", job_id, "/output/chr_1.zip")
  if(!file.exists(job_id_path)) next
  job_txt <- readLines(paste0("/pub59/iasiimwe/imputationserver/workspace/", job_id, "/logs/job.txt"))
  time_tb2[[i]] <- as.numeric(ymd_hms(gsub(" Cleanup.*", "", job_txt[length(job_txt)])) - ymd_hms(gsub(" Setup.*", "", job_txt[1])), "secs")
}
imputation <- tibble(time_s = time_tb2) %>%
  mutate(dataset = row_number(),
         Step = "Imputation on server")

# Minimac4 results extraction
results_extraction <- tibble(dataset = 1:n_datasets) %>%
  left_join(read_csv(paste0("michigan/time_", mechanism, "_", missing_percentage, "_", ref_population, ".csv"), show_col_types = FALSE)) %>%
  mutate(time_s = as.numeric(ymd_hms(end) - ymd_hms(start), "secs")) %>%
  select(time_s) %>%
  mutate(dataset = row_number(),
         Step = "Extracting imputed results")

# Info score extraction
info_score <- tibble(dataset = 1:n_datasets) %>%
  left_join(read_csv(paste0("michigan/info_scores_time_", ref_population, "_", mechanism, "_", missing_percentage, ".csv"), show_col_types = FALSE)) %>%
  mutate(time_s = as.numeric(ymd_hms(end) - ymd_hms(start), "secs")) %>%
  select(time_s) %>%
  mutate(dataset = row_number(),
         Step = "Extracting info scores")

results_tb <- monolix_fitting %>%
  bind_rows(minimac) %>%
  bind_rows(imputation) %>%
  bind_rows(results_extraction) %>%
  bind_rows(info_score)
write_csv(results_tb, "Michigan_detailed_time.csv")


# Generate plots
# --------------
for (effect_scenario in c("low", "high")) {
  results <- read_csv(paste0("performance_metrics_", effect_scenario, ".csv"), show_col_types = FALSE) %>%
    filter(!(Method == "michigan_afr_ms" & Percentage == 10))
  results_time <- read_csv(paste0("performance_time_", effect_scenario, ".csv"), show_col_types = FALSE) %>%
    filter(!(Method == "michigan_afr_ms" & Percentage == 10))
  
  # Plots
  # -----
  label_fn <- function(x) {
    x %>%
      gsub("pmm", "Predictive\nmean matching", .) %>%
      gsub("midastouch", "Weighted predictive\nmean matching", .) %>%
      gsub("sample", "Random sample from\nobserved values", .) %>%
      gsub("cart", "Classification and\nregression trees", .) %>%
      gsub("rf", "Random forest", .) %>%
      gsub("^mean$", "Unconditional\nmean", .) %>%
      gsub("norm", "Bayesian linear\nregression", .) %>%
      gsub("ri", "Random indicator for\nnonignorable data", .) %>%
      gsub("michigan_single", "Michigan cloned\nserver (single)", .) %>%
      gsub("michigan_ms", "Michigan cloned server\n(multiple, SNP exclusion)", .) %>%
      gsub("michigan_5dat", "Michigan cloned server\n(multiple, 25 datasets)", .) %>%
      gsub("michigan_50dat", "Michigan cloned server\n(multiple, 50 datasets)", .) %>%
      gsub("michigan", "Michigan cloned\nserver (multiple)", .) %>%
      gsub("CCA", "Complete Case\nAnalysis", .)
  } 
  
  results_success <- results %>%
    filter(Metric == "Success") %>%
    select(-Parameter, -Metric, -Label) %>%
    distinct() %>%
    mutate(Method = gsub("_afr", "", Method),
           Label = factor(Method, levels = unique(c("michigan_single", "CCA", "Mode", "michigan_ms", "michigan", Method))),
           Mechanism = gsub("_post", " (post-processing for MICE methods)", Mechanism),
           Mechanism = factor(Mechanism, levels = c("MCAR", "MAR", "MNAR", "MNAR (post-processing for MICE methods)"))) %>%
    filter(!Label %in% c("michigan_ms", "michigan_5dat", "michigan_50dat"))
  
  label_df <- results_success %>%
    group_by(Method, Mechanism, Label) %>%
    summarize(Value = max(Value), .groups = "drop")
  
  ggplot(data = filter(results_success, Percentage == 5), aes(x = Label, y = Value)) + 
    geom_bar(stat = "identity", alpha = 0.9, fill = brewer.pal(7,"Greens")[7]) +
    facet_wrap(~Mechanism, ncol = 1) +
    # 10%
    geom_bar(data = filter(results_success, Percentage == 10), 
             stat = "identity", width = 0.7, alpha = 0.7, fill = brewer.pal(7,"Greens")[4]) +
    # 20%
    geom_bar(data = filter(results_success, Percentage == 20), 
             stat = "identity", width = 0.5, alpha = 0.7, fill = brewer.pal(7,"Reds")[4]) +
    # 50%
    geom_bar(data = filter(results_success, Percentage == 50), 
             stat = "identity", width = 0.3, alpha = 0.7, fill = brewer.pal(7,"Reds")[7]) +
    # 100%
    geom_bar(data = filter(results_success, Label == "michigan_single"), 
             stat = "identity", alpha = 0.9, fill = "red") +
    # Add labels for 5%
    geom_text(data = label_df, 
              aes(label = ifelse(Value < 100, round(Value, 0), ""), 
                  vjust = ifelse(Value < 100, -0.2, 1.2))) +
    geom_text(data =label_df, 
              aes(label = ifelse(Value < 100, "", round(Value, 0)), 
                  vjust = ifelse(Value < 100, -0.2, 1.2)), colour = "white") +
    scale_x_discrete(labels = label_fn) +
    labs(title = NULL,
         x = "Imputation Method",
         y = "Number of Completed Runs", fill = "Percentage") +
    theme_bw() +
    theme(strip.text = element_text(size = 12, face = "bold"),
          axis.text = element_text(size = 10, colour = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.ticks.x = element_blank(), 
          axis.title = element_text(size = 12, face = "bold")) + 
    geom_vline(xintercept = 3.5, linetype = "dashed", color = "black", linewidth = 1) # Line to settle single and multiple imputation
  x <- 0.8
  if (file.exists(paste0("Fig_Success_", effect_scenario, ".png"))) file.remove(paste0("Fig_Success_", effect_scenario, ".png"))
  ggsave(paste0("Fig_Success_", effect_scenario, ".png"), height = 9 * x, width = 16 * x)
  
  
  # Create a data frame for the legend
  legend_data <- tibble(Percentage = c("5%", "10%", "20%", "50%", "100%"),
                        x = 1:5,  # Positions of the columns
                        y = 1,    # Single row
                        color = c(brewer.pal(7, "Greens")[7], 
                                  brewer.pal(7, "Greens")[4], 
                                  brewer.pal(7, "Reds")[4], 
                                  brewer.pal(7, "Reds")[7],
                                  "red")
  )
  
  # Create the heatmap-style legend
  ggplot(legend_data, aes(x = x, y = y, fill = color)) +
    geom_tile() +  # Draw tiles (heatmap blocks)
    geom_text(aes(label = Percentage), vjust = -1, size = 5, colour = "white") +  # Add labels above the tiles
    scale_fill_identity() +  # Use the colors directly without mapping them to a scale
    theme_void() +  # Remove axes and grid lines
    theme(legend.position = "none")  # Remove default legends
  x <- 0.2
  if (file.exists("Fig_Success_legend.png")) file.remove("Fig_Success_legend.png")
  ggsave("Fig_Success_legend.png", height = 9 * x, width = 16 * x, bg = "white")
  
  
  # Performance metrics
  # eMLAR (Bias) and eMALAR (Precision) unbiased - and conclusions similar to the rest, so shown
  if (!dir.exists(paste0(effect_scenario, "_effect_plots"))) dir.create(paste0(effect_scenario, "_effect_plots"))
  
  label_fn <- function(x) {
    x %>%
      gsub("pmm", "Predictive\nmean matching", .) %>%
      gsub("midastouch", "Weighted predictive\nmean matching", .) %>%
      gsub("sample", "Random sample from\nobserved values", .) %>%
      gsub("cart", "Classification and\nregression trees", .) %>%
      gsub("rf", "Random forest", .) %>%
      gsub("^mean$", "Unconditional\nmean", .) %>%
      gsub("norm", "Bayesian linear\nregression", .) %>%
      gsub("ri", "Random indicator for\nnonignorable data", .) %>%
      gsub("michigan_single", "Michigan cloned\nserver (single)", .) %>%
      gsub("michigan_ms", "Michigan cloned server\n(multiple, SNP exclusion)", .) %>%
      gsub("michigan_5dat", "Michigan cloned server\n(multiple, 2 datasets)", .) %>%
      gsub("michigan_50dat", "Michigan cloned server\n(multiple, 50 datasets)", .) %>%
      gsub("michigan", "Michigan cloned\nserver (multiple)", .) %>%
      gsub("CCA", "Complete Case\nAnalysis", .)
  } 
  
  performance_plot_fn <- function(dat, metric = "eMLAR", params = c("Cl THETA", "Cl ETA"),
                                  subset_methods = NULL) {
    dat_plot <- dat %>% 
      filter(!is.na(Value)) %>%
      filter(Metric %in% metric, Parameter %in% params) %>%
      filter(!str_detect(Parameter, "RSE")) %>%
      mutate(Parameter = gsub("ka THETA", "Absorption rate constant (/hr)", Parameter),
             Parameter = gsub("ka ETA", "Between-subject variability, absorption rate constant (%)", Parameter),
             Parameter = gsub("Cl THETA", "Apparent clearance (L/hr)", Parameter),
             Parameter = gsub("Cl ETA", "Between-subject variability, Clearance (%)", Parameter),
             Parameter = gsub("^b$", "Proportional error (%)", Parameter),
             Parameter = factor(Parameter, levels = unique(Parameter)),
             Method = gsub("_afr", "", Method),
             Mechanism = gsub("_post", " (post-processing for MICE methods)", Mechanism),
             Mechanism = factor(Mechanism, levels = c("MCAR", "MAR", "MNAR")),
             Fill = if_else(str_detect(Method, "michigan_single|michigan_ms"), "100", as.character(Percentage)),
             Fill = factor(Fill, levels = c("5", "10", "20", "50", "100"))
      )
    method_levels <- unique(c("michigan_single", "CCA", "Mode", "michigan_ms", "michigan", dat_plot$Method))
    if (!is.null(subset_methods)) method_levels <- method_levels[method_levels %in% subset_methods]
    if (!is.null(subset_methods)) dat_plot <- filter(dat_plot, Method %in% subset_methods) 
    # %>% mutate(Fill = factor(as.character(Percentage), levels = c("5", "10", "20", "50")))
    dat_plot <- mutate(dat_plot, Method = factor(Method, levels = method_levels))
    
    # Define custom facet labels
    custom_labeller <- function(x) {
      x[duplicated(x)] <- ""  # Remove repeated labels
      return(x)
    }
    
    dat_plot <- dat_plot %>%
      ggplot(aes(x = Method, y = Value, fill = Fill)) +
      geom_bar(stat = "identity", position = position_dodge2(preserve = "single", padding = 0), width = 0.7) +
      facet_nested_wrap(vars(Mechanism, Parameter), nrow = 4, scales = "free_y",
                        labeller = labeller(Parameter = custom_labeller)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
      scale_fill_manual(values = c("5" = brewer.pal(7,"Greens")[7], "10" = brewer.pal(7,"Greens")[4], 
                                   "20" = brewer.pal(7,"Reds")[4], "50" = brewer.pal(7,"Reds")[7], "100" = "red")) +
      scale_x_discrete(labels = label_fn) +
      labs(title = NULL,
           x = "Imputation Method",
           y = "Percentage", fill = "Percentage") +
      theme_bw() + 
      theme(legend.position = "bottom", 
            legend.background = element_rect(fill = "white", colour = "black"), # Optional: Add a white background to the legend
            strip.text = element_text(size = 12, face = "bold", margin = margin(0.05, 0, 0.05, 0, "cm")),
            axis.text = element_text(size = 11, colour = "black"),
            strip.background = element_blank(),  # Removes strip background for better blending
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title = element_text(size = 12, face = "bold"),
            plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) 
    
    if (metric == "CR") dat_plot <- dat_plot + 
      geom_hline(yintercept = 90, linetype = "dashed", color = "black") +
      coord_cartesian(ylim = c(0, 100))
    if (metric == "AW") dat_plot <- dat_plot + ylab("Confidence Interval Width") 
    return(dat_plot)
  }
  
  
  metrics <- c("eMLAR", "MRPE", "rMPE", "eMALAR", "RMSRE", "MAPE", "AW", "CR")
  # params_list <- list(c("Cl THETA", "Cl ETA"), c("ka THETA", "ka ETA"), c("b"))
  params_list <- list(c("Cl THETA"), c("Cl ETA"), c("ka THETA"), c("ka ETA"), c("b"))
  for (metric in metrics) {
    for (params in params_list) {
      params <- unlist(params)
      plot <- performance_plot_fn(results, metric, params) +
        geom_text(aes(label = round(Value, 0)),
                  position = position_dodge2(preserve = "single", padding = 0, width = 0.9),
                  vjust = -0.2, size = 3.5) +
        geom_vline(xintercept = 3.5, linetype = "dashed", color = "black", linewidth = 1) 
      # + ggtitle(gsub(" ", "_", paste0(params, collapse = "_")))
      file_name <- paste0(effect_scenario, "_effect_plots/", metric, "_", effect_scenario)
      if (file.exists(file_name)) file.remove(file_name)
      x <- 1
      # ggsave(file_name, height = 9 * x, width = 16 * x)
      graph2ppt(plot, file = file_name, append = TRUE, width = 16 * x, height = 9 * x)
    }
  }
  
  
  # Time
  # -------
  # Michigan detailed time
  michigan <- read_csv("Michigan_detailed_time.csv", show_col_types = FALSE)
  datasets_with_na <- michigan %>%
    filter(is.na(time_s)) %>%
    pull(dataset) %>%
    unique()
  length(datasets_with_na) # 44, exclude these as imputation wasn't successful
  michigan <- michigan %>%
    filter(!dataset %in% datasets_with_na) %>%
    mutate(Step = factor(Step, levels = c("Preparation for imputation", "Imputation on server", "Extracting imputed results",
                                          "Extracting info scores", "Model fitting (Monolix)")),
           time_m = time_s/60,
           `Step Similarity (Single vs. Multiple Imputation)` = ifelse(Step == "Model fitting (Monolix)", "Different", "Similar"))
  
  # Plot
  ggplot(michigan, aes(x = Step, y = time_m, colour = `Step Similarity (Single vs. Multiple Imputation)`)) +
    geom_boxplot() +
    coord_flip() +  # Flips the axes for better readability
    labs(title = "Distribution of Time Spent on Each Step",
         x = "Step",
         y = "Time (minutes)") +
    theme_bw() +
    theme(legend.position = "bottom", 
          axis.text = element_text(size = 10, colour = "black"),
          axis.title = element_text(size = 12, face = "bold"),
          plot.title = element_text(hjust = 0, size = 14, face = "bold")) 
  x <- 0.45
  ggsave("Michigan_detailed_time.png", height = 9 * x, width = 16 * x, bg = "white")
  
  methods_a <- c("CCA", "Mode", "michigan_afr_single", "Reference") # Don't show multiple imputation because it is similar in time (see results above)
  methods_b <- c("pmm", "midastouch", "sample", "cart")
  methods_c <- c("rf", "mean", "norm", "ri")
  label_fn <- function(x) {
    x <- if_else(str_detect(x, "Ref|mic|10"), gsub("5|10", "", x), " ") %>%
      gsub("michigan_afr_single.*", "Michigan (single)", .) %>%
      gsub("michigan_afr.*", "Michigan (multiple)", .) %>%
      gsub("midastouch*", "Weighted\npmm", .) 
    return(x)
  }
  time_plot_fn <- function(dat, y_break = c(10, 500), to_break = TRUE) {
    dat <- dat %>%
      bind_rows(mutate(distinct(filter(dat, Method == "Reference", Percentage == 5)), Mechanism = "MAR")) %>%
      bind_rows(mutate(distinct(filter(dat, Method == "Reference", Percentage == 5)), Mechanism = "MNAR")) %>%
      filter(!(str_detect(Method, "michigan") & Percentage %in% c(20, 50))) %>%
      filter(!is.na(time_s)) %>%
      mutate(Method2 = paste0(Method, Percentage),
             Method2 = factor(Method2, levels = unique(c("Reference5", "michigan_afr_single5", Method2))),
             Method = gsub("michigan_", "Michigan ", Method),
             Method = gsub("afr_single.*", "(single)", Method),
             Method = gsub("afr.*", "(multiple)", Method),
             Method = gsub("Reference.*", "Original", Method),
             Method = factor(Method, levels = unique(c("Original", "Michigan (single)", Method))),
             Mechanism = gsub("_post", " (post)", Mechanism),
             Mechanism = factor(Mechanism, levels = c("MCAR", "MAR", "MNAR")),
             Fill = if_else(!str_detect(Method, "Michigan|Original"), as.character(Percentage), Method),
             Fill = gsub("Original", "N/A", Fill),
             Fill = gsub("Michigan \\(single\\)", "All", Fill),
             Fill = factor(Fill, levels = c("N/A", "All", "5", "10", "20", "50")))
    median_labels <- dat %>% # Calculate medians for each Method, Percentage and Mechanism combination
      group_by(Method2, Mechanism) %>%
      summarise(median_time = median(time_s), .groups = "drop") %>%
      mutate(Label = if_else(median_time < 60, paste0(round(median_time, 1), " s"), paste0(round(median_time/60, 1), " m")), # Convert seconds to minutes
             Label = if_else(str_detect(Label, "\\."), Label, gsub(" ", ".0 ", Label)),
             Label = gsub(" ", "", Label),
             hjust_dynamic = ifelse(median_time > 40, 1.5, -2.1)) # To add median labels dynamically
    dat <- dat %>%
      ggplot(aes(x = Method2, y = time_s)) +
      geom_boxplot(aes(colour = Fill)) +  
      geom_text(data = median_labels, 
                aes(x = Method2, y = median_time, label = Label, hjust = hjust_dynamic), size = 3) +  
      facet_wrap(~ Mechanism, ncol = 1, strip.position="right")
    if(to_break) dat <- dat + scale_y_break(y_break, scales = c(1, 1))
    dat +
      scale_y_log10() +  # Apply log scale to y-axis
      scale_x_discrete(labels = label_fn) +
      scale_colour_manual(values = c(brewer.pal(3,"Greens")[3], brewer.pal(3,"Blues")[3], brewer.pal(7,"Reds")[c(2, 3, 5, 7)])) +
      labs(title = NULL,
           x = "Method",
           y = "Time (Seconds)",
           colour = "Percentage") +
      theme_bw() + 
      theme(legend.position = "none",
            strip.text = element_text(size = 12, face = "bold"),
            axis.text = element_text(size = 10, colour = "black"),
            axis.text.x = element_text(hjust = 1),
            axis.ticks.y = element_blank(), 
            axis.title = element_text(size = 12, face = "bold"),
            plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
      coord_flip() 
  }
  # time_plot_fn(results_time, Methods = c("CCA", "Mode"), y_break = c(15, 90), to_break = FALSE)
  a <- time_plot_fn(results_time %>% filter(Method %in% methods_a), y_break = c(10, 1000)) + ylab("   ")
  b <- time_plot_fn(results_time %>% filter(Method %in% methods_b), to_break = FALSE) + xlab("   ") +
    theme(legend.position = "top")
  c <- time_plot_fn(results_time %>% filter(Method %in% methods_c), to_break = FALSE) + xlab("   ") + ylab("   ") 
  a
  plot_name <- paste0("Time1_", effect_scenario, ".png")
  if (file.exists(plot_name)) file.remove(plot_name)
  x <- 0.9
  ggsave(plot_name, height = 9 * x, width = 8 * x)
  b
  plot_name <- paste0("Time2_", effect_scenario, ".png")
  if (file.exists(plot_name)) file.remove(plot_name)
  x <- 0.8
  ggsave(plot_name, height = 9 * x, width = 10 * x)
  c
  plot_name <- paste0("Time3_", effect_scenario, ".png")
  if (file.exists(plot_name)) file.remove(plot_name)
  x <- 0.75
  ggsave(plot_name, height = 9 * x, width = 9 * x)
}



