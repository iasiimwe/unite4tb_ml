# Load required packages
# ---------------------
library(lixoftConnectors)
initializeLixoftConnectors(software = "monolix", path = "/pub59/iasiimwe/Lixoft/MonolixSuite2024R1/", force = TRUE)
library(tidyverse)
library(data.table)

select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate

# Add folders if they don't exist
if (!dir.exists("missingness")) dir.create("missingness")

# Dominant coding function
dominant_coding_fn <- function(x) {
  x <- gsub(" ", "", x) 
  y <- table(x)
  y <- names(sort(-y)) # Sort in descending order, assume mutant type is the least abundant
  if (length(y) == 3) x <- gsub(y[[3]], "1", x) # Give mutant homozygotes "1"
  if (length(y) > 1) x <- gsub(y[[2]], "1", x) # Give heterozygotes "1" 
  x <- gsub(y[[1]], "0", x) # Give homozygotes 0
  return(as.numeric(x))
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

message(paste0("Starting complete data analysis!"))

# Results storage matrix
params_for_estimation <- c("ka_theta", "ka_theta_rse", "ka_eta", "ka_eta_rse", 
                           "Cl_theta", "Cl_theta_rse", "Cl_eta", "Cl_eta_rse", "b", "b_rse")
results_matrix <- matrix(0, nrow = n_datasets, ncol = length(params_for_estimation),
                         dimnames = list(paste0("dataset_", 1:n_datasets), params_for_estimation))
ml_methods <- "original_data"
time_list <- rep(list(tibble(end = vector("character", n_datasets),
                             start = vector("character", n_datasets))), 
                 length(ml_methods))
k <- 1 # Only one method

for (j in 1:n_datasets) {
  dat <- read_csv(paste0("sim_DV/sim_data", j, "_", sim_date, ".csv"), show_col_types = FALSE) %>%
    select(-NTIME, -ID2) 

  # Get time of code execution start
  start <- Sys.time()

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
  time_list[[k]]$start[j] <- paste("T", as.character(ymd_hms(start)))
  time_list[[k]]$end[j] <- paste("T", as.character(ymd_hms(end)))
  
  message(paste0("##################\n##################\n##################\n", 
                 round(j * 100/n_datasets, 2), 
                 "% complete\n##################\n##################\n##################"))
}
write.csv(results_matrix, paste0("missingness/complete.csv"), row.names = FALSE, quote = FALSE)
write.csv(time_list[[k]], paste0("missingness/", ml_methods[[k]], "_time.csv"), row.names = FALSE)

# Obtain plots
# ------------
# Reshape the data for plotting
dat <- read_csv(paste0("missingness/complete.csv"), show_col_types = FALSE)
colnames(dat) <- gsub("b", "b_b", gsub("_rse", "rse", colnames(dat)))

dat_long <- dat %>%
  pivot_longer(
    cols = everything(),
    names_to = c("parameter", "metric"),
    names_sep = "_",
    values_to = "value"
  ) %>%
  mutate(metric = toupper(metric),
         metric = gsub("ETA", "ETAs", metric),
         metric = gsub("RSE", " (RSEs)", metric),
         metric = gsub("B", "b", metric),
         metric = factor(metric, levels = unique(metric)),
         parameter = factor(parameter, levels = unique(parameter))) %>% 
  filter(value < 100) # 9 eta RSEs greater than 100 not shown


# Create the boxplots with facets for metrics (theta, eta, shrinkage)
ggplot(dat_long, aes(x = parameter, y = value, colour = parameter)) +
  geom_boxplot() +
  facet_wrap(~metric, scales = "free", nrow = 3) +
  theme_bw() +
  labs(
    title = "Boxplots of Estimates for Parameters Ka, Cl, and proportional error (b)",
    x = "Parameter",
    y = "Value"
  ) +
  theme(legend.position = "none",
        strip.text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10, colour = "black"),
        axis.text.x = element_text(angle = 0, hjust = 1, size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0, vjust = 0, size = 14, face = "bold"))
x <- 0.7
ggsave("true_estimates.png", height = 12 * x, width = 9 * x)







