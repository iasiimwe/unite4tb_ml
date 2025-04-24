# Load required packages
# ---------------------
library(tidyverse)
library(ggh4x)
library(ggrepel)
library(lubridate)

# Relevant functions
ml_mtd_tb <- function(method_tb, method, col_to_use = "Correlated_covar1", nsims = 100) {
  
  # Add missing Simulations
  included_sims <- unique(method_tb$Simulation)
  missing_sims <- c(1:nsims)[!c(1:nsims) %in% included_sims]
  if (length(missing_sims > 0)) {
    missing_sims <- tibble(Method = method,
                           Simulation = missing_sims,
                           Covariates = "",
                           Correlated_covar1 = "",
                           Correlated_covar2 = "",
                           Importance = "")
    method_tb <- method_tb %>%
      bind_rows(missing_sims) %>%
      arrange(Simulation, desc(Importance))
  }
  
  method_tb$Correlated_covar <- method_tb[[col_to_use]]
  
  summary_tb <- method_tb %>%
    group_by(Simulation) %>%
    summarise(top_m = paste0(Covariates[1:n_true_snps], collapse = ", "),
              top_m = gsub(", NA|NA, |NA", "", top_m),
              Correlated_covar = gsub("^, |, $| $", "", paste0(unique(unlist(str_split(Correlated_covar[1:n_true_snps], ", "))), collapse = ", ")),
              Correlated_covar = gsub(", ,", ",", Correlated_covar),
              Correlated_covar = gsub(", NA|NA, |NA", "", Correlated_covar)) %>%
    rowwise() %>%
    mutate(Correlated_covar = paste0(unique(c(unlist(str_split(top_m, ", ")), unlist(str_split(Correlated_covar, ", ")))), collapse = ", "),
           n_correlated_covar = str_count(Correlated_covar, ",") + 1,
           n_covar = str_count(top_m, ",") + 1, 
           TruePositive = sum(str_detect(top_m, true_covar)),
           TruePositivecor = sum(str_detect(paste0(top_m, Correlated_covar), true_covar)),
           percent_correct = sum(str_detect(top_m, true_covar))*100/length(true_covar),
           percent_correctcor = sum(str_detect(paste0(top_m, Correlated_covar), true_covar))*100/length(true_covar),
           detected_true_covar = paste0(true_covar[str_detect(top_m, true_covar)], collapse = ", "),
           detected_true_covarcor = paste0(true_covar[str_detect(paste0(top_m, Correlated_covar), true_covar)], collapse = ", "),
           detected_false_covar = paste0(setdiff(unlist(str_split(top_m, ", ")), 
                                                 unlist(str_split(detected_true_covar, ", "))), 
                                         collapse = ", ")) %>%
    ungroup()
  
  true_covar_table <- tibble(Covariate = true_covar,
                             times = 0)
  for (l in seq_along(true_covar)) true_covar_table$times[[l]]  <- sum(str_detect(summary_tb$detected_true_covar, true_covar[[l]]))
  
  cor_covar_table <- tibble(Covariate = true_covar,
                            times = 0)
  for (l in seq_along(true_covar)) cor_covar_table$times[[l]]  <- sum(str_detect(summary_tb$Correlated_covar, true_covar[[l]]))
  
  # Top-M selection
  top_m <- summary_tb  %>%
    summarize(Method = method,
              `Average % correct` = round(mean(percent_correct), n_round),
              `Average % correctcor` = round(mean(percent_correctcor), n_round),
              n_correlated_covar_median = round(median(n_correlated_covar), n_round),
              n_correlated_covar_min = round(min(n_correlated_covar), n_round),
              n_correlated_covar_max = round(max(n_correlated_covar), n_round),
              n_covar_median = round(median(n_covar), n_round),
              n_covar_min = round(min(n_covar), n_round),
              n_covar_max = round(max(n_covar), n_round),
              `Perfect selection %` = round(sum(percent_correct == 100)*100/nsims, n_round),
              `Perfect selectioncor %` = round(sum(percent_correctcor == 100)*100/nsims, n_round),
              `True covariate % selection` = paste(true_covar_table$Covariate, round(true_covar_table$times*100/nsims, n_round), sep = ", ", collapse = "; "),
              `Correlated covariate % selection` = paste(cor_covar_table$Covariate, round(cor_covar_table$times*100/nsims, n_round), sep = ", ", collapse = "; "))
  
  f1 <- summary_tb  %>%
    mutate(FalsePositive = n_true_snps - TruePositive,
           FalsePositivecor = n_true_snps - TruePositivecor,
           FalseNegative = length(true_covar) - TruePositive,
           FalseNegativecor = length(true_covar) - TruePositivecor,
           Recall = TruePositive/ (TruePositive + FalseNegative),
           Recallcor = TruePositivecor/ (TruePositivecor + FalseNegativecor),
           Precision = TruePositive/ (TruePositive + FalsePositive),
           Precisioncor = TruePositivecor/ (TruePositivecor + FalsePositivecor),
           `F1 score` = 2 * (Recall * Precision)/(Recall + Precision),
           `F1 score cor` = 2 * (Recallcor * Precisioncor)/(Recallcor + Precisioncor),
           `F1 score` = ifelse(is.na(`F1 score`) == TRUE, 0, `F1 score`),
           `F1 score cor` = ifelse(is.na(`F1 score cor`) == TRUE, 0, `F1 score cor`)) %>%
    summarise(`F1 score mean` = mean(`F1 score`),
              `F1 score sd` = sd(`F1 score`),
              `F1 score cor mean` = mean(`F1 score cor`),
              `F1 score cor sd` = sd(`F1 score cor`)) # Will need to adapt this to multiple scenarios and plotting ROC curves
  
  return(bind_cols(top_m, f1))
}

summary_correlated <- function(path) {
  read_csv(path, show_col_types = FALSE) %>%
    select(Method, n_correlated_covar_median, n_correlated_covar_min, n_correlated_covar_max,
           `Perfect selectioncor %`, `F1 score cor mean`, `F1 score cor sd`, Scenario, NSNPs) %>%
    arrange(Scenario, NSNPs, desc(`F1 score cor mean`), desc(`Perfect selectioncor %`)) %>%
    mutate(n = paste0(n_correlated_covar_median, " (", n_correlated_covar_min, " to ", n_correlated_covar_max, ")"),
           F1 = paste0(round(`F1 score cor mean`, n_round), " (", round(`F1 score cor sd`, n_round), ")"),
           perfect = `Perfect selectioncor %`) %>%
    select(Method, F1, perfect, n, Scenario, NSNPs) %>%
    return()
}

summary_not_correlated <- function(path) {
  read_csv(path, show_col_types = FALSE) %>%
    select(Method, n_covar_median, n_covar_min, n_covar_max, 
           `Perfect selection %`, `F1 score mean`, `F1 score sd`, Scenario, NSNPs) %>%
    arrange(Scenario, NSNPs, desc(`F1 score mean`), desc(`Perfect selection %`)) %>%
    mutate(n = paste0(n_covar_median, " (", n_covar_min, " to ", n_covar_max, ")"),
           F1 = paste0(round(`F1 score mean`, n_round), " (", round(`F1 score sd`, n_round), ")"),
           perfect = `Perfect selection %`) %>%
    select(Method, F1, perfect, n, Scenario, NSNPs) %>%
    return()
}

# General settings
n_datasets <- 100
scenarios <- c("low", "high")
true_covar <- paste0("SNP", 1:9)
n_true_snps <- length(true_covar)
n_round <- 3 
n_snp <- "10_3"
ml_methods <- c("glmnet", "rf", "gwas", "gwas_glmnet", "gwas_rf") 

performance_tb1 <- tibble()
performance_tb2 <- tibble()
for (scenario in scenarios) {
  for (ml_method in ml_methods) {
    message(paste0("Scenario: ", scenario, 
                   "\nMethod: ", ml_method))
    method_tb_path <- paste0(scenario, "/covariates/", n_snp, "_", ml_method, "_sim_results.csv")
    method_tb <- read_csv(method_tb_path, show_col_types = FALSE) %>%
      mutate(Covariates = gsub("`", "", Covariates),
             Correlated_covar1 = gsub("`", "", Correlated_covar1),
             Correlated_covar2 = gsub("`", "", Correlated_covar2),
             Correlated_covar1 = ifelse(is.na(Correlated_covar1), "", Correlated_covar1),
             Correlated_covar2 = ifelse(is.na(Correlated_covar2), "", Correlated_covar2))
    
    performance_tb1 <- bind_rows(performance_tb1, 
                                 ml_mtd_tb(method_tb, ml_method, col_to_use = "Correlated_covar1") %>%
                                   mutate(Scenario = scenario, NSNPs = n_snp))
    performance_tb2 <- bind_rows(performance_tb2, 
                                 ml_mtd_tb(method_tb, ml_method, col_to_use = "Correlated_covar2") %>%
                                   mutate(Scenario = scenario, NSNPs = n_snp))
  }
}
write.csv(performance_tb1, "Eta_sim_performance1_1000.csv", row.names = FALSE)
write.csv(performance_tb2, "Eta_sim_performance2_1000.csv", row.names = FALSE)

# Combined
# --------
# Not_correlated
summary_not_correlated("Eta_sim_performance1_1000.csv") %>%
  write.csv(., "eta_not_correlated_1000.csv", row.names = FALSE)

# Correlated
summary_correlated("Eta_sim_performance1_1000.csv") %>%
  write.csv(., "eta_correlated1_1000.csv", row.names = FALSE)
summary_correlated("Eta_sim_performance2_1000.csv") %>%
  write.csv(., "eta_correlated2_1000.csv", row.names = FALSE)

