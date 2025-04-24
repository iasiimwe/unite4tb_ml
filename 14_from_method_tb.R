# Load required packages
# ---------------------
library(tidyverse)
library(ggh4x)
library(ggrepel)
library(lubridate)

# Relevant functions
ml_mtd_tb <- function(method_tb, method, nsims = 100) {
  
  # Add missing Simulations
  included_sims <- unique(method_tb$Simulation)
  missing_sims <- c(1:nsims)[!c(1:nsims) %in% included_sims]
  if (length(missing_sims > 0)) {
    missing_sims <- tibble(Method = method,
                           Simulation = missing_sims,
                           Covariates = "",
                           Importance = "")
    method_tb <- method_tb %>%
      bind_rows(missing_sims) %>%
      arrange(Simulation, desc(Importance))
  }

  summary_tb <- method_tb %>%
    group_by(Simulation) %>%
    summarise(top_m = paste0(Covariates[1:n_true_snps], collapse = ", "),
              top_m = gsub(", NA|NA, |NA", "", top_m)) %>%
    rowwise() %>%
    mutate(n_covar = str_count(top_m, ",") + 1, 
           TruePositive = sum(str_detect(top_m, true_covar)),
           percent_correct = sum(str_detect(top_m, true_covar))*100/length(true_covar),
           detected_true_covar = paste0(true_covar[str_detect(top_m, true_covar)], collapse = ", "),
           detected_false_covar = paste0(setdiff(unlist(str_split(top_m, ", ")), 
                                                 unlist(str_split(detected_true_covar, ", "))), 
                                         collapse = ", ")) %>%
    ungroup()
  
  true_covar_table <- tibble(Covariate = true_covar,
                             times = 0)
  for (l in seq_along(true_covar)) true_covar_table$times[[l]]  <- sum(str_detect(summary_tb$detected_true_covar, true_covar[[l]]))
  
  # Top-M selection
  top_m <- summary_tb  %>%
    summarize(Method = method,
              `Average % correct` = round(mean(percent_correct), n_round),
              n_covar_median = round(median(n_covar), n_round),
              n_covar_min = round(min(n_covar), n_round),
              n_covar_max = round(max(n_covar), n_round),
              `Perfect selection %` = round(sum(percent_correct == 100)*100/nsims, n_round),
              `True covariate % selection` = paste(true_covar_table$Covariate, round(true_covar_table$times*100/nsims, n_round), sep = ", ", collapse = "; "))
  
  f1 <- summary_tb  %>%
    mutate(FalsePositive = n_true_snps - TruePositive,
           FalseNegative = length(true_covar) - TruePositive,
           Recall = TruePositive/ (TruePositive + FalseNegative),
           Precision = TruePositive/ (TruePositive + FalsePositive),
           `F1 score` = 2 * (Recall * Precision)/(Recall + Precision),
           `F1 score` = ifelse(is.na(`F1 score`) == TRUE, 0, `F1 score`)) %>%
    summarise(`F1 score mean` = mean(`F1 score`),
              `F1 score sd` = sd(`F1 score`)) 
  
  return(bind_cols(top_m, f1))
}

summary_not_correlated <- function(path) {
  read_csv(path, show_col_types = FALSE) %>%
    select(Method, n_covar_median, n_covar_min, n_covar_max, 
           `Perfect selection %`, `F1 score mean`, `F1 score sd`, Scenario, Pruned, NSNPs) %>%
    arrange(Scenario, Pruned, NSNPs, desc(`F1 score mean`), desc(`Perfect selection %`)) %>%
    mutate(n = paste0(n_covar_median, " (", n_covar_min, " to ", n_covar_max, ")"),
           F1 = paste0(round(`F1 score mean`, n_round), " (", round(`F1 score sd`, n_round), ")"),
           perfect = `Perfect selection %`) %>%
    select(Method, F1, perfect, n, Scenario, Pruned, NSNPs) %>%
    return()
}

# General settings
n_datasets <- 100
scenarios <- c("low", "high") 
true_covar <- paste0("SNP", 1:9)
n_true_snps <- length(true_covar)
n_round <- 3 
n_round_time <- 2
n_snps <- c("10_3", "10_4", "10_5", "10_6")
pruned_snps <- c("none", "pruned", "pca", "pca_o")
ml_methods <- c("glmnet", "rf", "gwas", "gwas_glmnet", "gwas_rf")
performance_tb1 <- tibble()
performance_tb2 <- tibble()
performance_tb_plot1 <- tibble()
performance_tb_plot2 <- tibble()
time_tb <- tibble()

for (scenario in scenarios) {
  for (pruned in pruned_snps) {
     for (n_snp in n_snps) {
      if (str_detect(pruned, "pca") & n_snp == "10_6") next
      for (ml_method in ml_methods) {
        if (str_detect(pruned, "pca") & str_detect(ml_method, "gwas")) next
        if (str_detect(pruned, "pruned") & !str_detect(ml_method, "gwas") & n_snp %in% c("10_6")) next
        if (str_detect(pruned, "none") & !str_detect(ml_method, "gwas") & n_snp %in% c("10_5", "10_6")) next
        
        message(paste0("Scenario: ", scenario, 
                       "\nPruned: ", pruned,
                       "\nN SNPs: ", n_snp,
                       "\nMethod: ", ml_method))
        
        method_tb_path <- paste0(scenario, "/covariates/", n_snp, "_", ml_method, "_sim_results.csv")
        if (pruned == "pruned") method_tb_path <- gsub("covariates/", "covariates/pruned_", method_tb_path)
        if (pruned == "pca") method_tb_path <- gsub("covariates/", "covariates/pca_", method_tb_path)
        if (pruned == "pca_o") method_tb_path <- gsub("covariates/", "covariates/pca_o_", method_tb_path)

        method_tb <- read_csv(method_tb_path, show_col_types = FALSE) %>%
          mutate(Covariates = gsub("`", "", Covariates))
        
        performance_tb1 <- bind_rows(performance_tb1, 
                                     ml_mtd_tb(method_tb, ml_method) %>%
                                       mutate(Scenario = scenario, Pruned = pruned, NSNPs = n_snp))

        # Time
        time_path <- paste0(scenario, "/covariates/", n_snp, "_", ml_method, "_time.csv")
        if (pruned == "pruned") time_path <- gsub("covariates/", "covariates/pruned_", time_path)
        if (pruned == "pca") time_path <- gsub("covariates/", "covariates/pca_", time_path)
        if (pruned == "pca_o") time_path <- gsub("covariates/", "covariates/pca_o_", time_path)
        
        time_s <- read_csv(time_path, show_col_types = FALSE) %>%
          mutate(end = ymd_hms(gsub("T ", "", end)),
                 start = ymd_hms(gsub("T ", "", start)),
                 diff = as.numeric(end-start, "secs")) %>%
          pull(diff)
        
        if (str_detect(ml_method, "gwas")) {
          gwas_time_path <- paste0(scenario, "/gwas/", n_snp, "_time_gwas.csv")
          if (pruned %in% c("pruned")) gwas_time_path <- gsub("gwas/", "gwas/pruned_", gwas_time_path)

          time_s2 <- read_csv(gwas_time_path, show_col_types = FALSE) %>%
            mutate(end = ymd_hms(gsub("T ", "", end)),
                   start = ymd_hms(gsub("T ", "", start)),
                   diff = as.numeric(end-start, "secs")) %>%
            pull(diff)
          time_s <- time_s + time_s2
        }
        
        time_tb <- bind_rows(time_tb, 
                             tibble(Method = ml_method,
                                    mean = mean(time_s), 
                                    sd = sd(time_s), 
                                    Scenario = scenario, 
                                    Pruned = pruned, 
                                    NSNPs = n_snp))
      }
    }
  }
}
write.csv(performance_tb1, "Eta_sim_performance1.csv", row.names = FALSE)
write.csv(performance_tb2, "Eta_sim_performance2.csv", row.names = FALSE)
write.csv(time_tb, "Time.csv", row.names = FALSE)

# Combined
# --------
# Not_correlated
summary_not_correlated("Eta_sim_performance1.csv") %>%
  write.csv(., "eta_not_correlated.csv", row.names = FALSE)

# Table 1
summary_not_correlated("Eta_sim_performance1.csv") %>%
  mutate(Scenario = factor(Scenario, levels = c("low", "high")),
         Pruned = factor(Pruned, levels = c("none", "pruned", "pca", "pca_o")),
         NSNPs = factor(NSNPs, levels = unique(NSNPs)),
         F12 = as.numeric(gsub(" .*", "", F1)),
         Method = case_when(Method == "gwas" ~ "GWAS",
                            Method == "rf" ~ "Random Forest",
                            Method == "glmnet" ~ "Penalized regression",
                            Method == "gwas_glmnet" ~ "GWAS + Penalized regression",
                            Method == "gwas_rf" ~ "GWAS + Random Forest")) %>%
  arrange(Scenario, Pruned, NSNPs, F12) %>%
  select(-n, -F12, -perfect) %>%
  pivot_wider(id_cols = c(Method, NSNPs),  # columns to keep as-is
              names_from = c(Scenario, Pruned),  # each combination will become a separate column
              values_from = F1) %>%
  write.csv("Table1.csv", row.names = FALSE)


