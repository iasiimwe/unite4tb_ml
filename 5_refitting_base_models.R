# Load required packages
# ---------------------
library(lixoftConnectors)
initializeLixoftConnectors(software = "monolix", path = "C:/Program Files/Lixoft/MonolixSuite2024R1/", force = TRUE)
library(tidyverse)
library(data.table)

select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate

# Add folders if they don't exist
if (!dir.exists("base_model_refit")) dir.create("base_model_refit")

## Random effects (considering log distribution)
omega_squared_fn <- function(percent_CV) log((percent_CV/100)^2 + 1)

# # Without re-inflation
omega_ka <- omega_squared_fn(68.5)
omega_Cl <- omega_squared_fn(17.5)
omega_V <- omega_squared_fn(18.2)
# gamma_V <- omega_squared_fn(11.7)

sim_date <-"23_12_2024"
n_datasets <- 100
effect_size <- 0.5

# Results storage matrix
params_for_estimation <- c("ka_theta", "ka_theta_rse", "ka_eta", "ka_eta_rse", "ka_eta_shr", 
                           "Cl_theta", "Cl_theta_rse", "Cl_eta", "Cl_eta_rse", "Cl_eta_shr",
                           "V_theta", "V_theta_rse", "V_eta", "V_eta_rse", "V_eta_shr")
results_matrix <- matrix(0, nrow = n_datasets, ncol = length(params_for_estimation),
                         dimnames = list(paste0("dataset_", 1:n_datasets), params_for_estimation))

for (j in 1:n_datasets) {
  dat <- read_csv(paste0("sim_DV/sim_data", j, "_", sim_date, ".csv"), show_col_types = FALSE) %>%
    select(!c(contains("SNP"), NTIME, ID2))
  write.csv(dat, "temp.csv", row.names = FALSE)

  initializeLixoftConnectors(software = "monolix", path = "C:/Program Files/Lixoft/MonolixSuite2024R1/", force = TRUE)
  newProject(data = list(dataFile = 'temp.csv',
                         headerTypes = c("id", "time", "evid", "occ", "amount", "observation", 
                                         "catcov", "regressor")),
             modelFile = 'tb_base_vinnard.txt')
  
  # Change error model, distribution, variability and correlations
  setErrorModel(DV = "combined1")
  setIndividualParameterDistribution(ka = "logNormal", TVV = "logNormal", TVCl = "logNormal")
  setIndividualParameterVariability(id = list (ka = TRUE, TVV = TRUE, TVCl = TRUE),
                                    OCC = list (ka = FALSE, TVV = FALSE, TVCl = FALSE))
  # getIndividualParameterModel() to confirm changes

  # getPopulationParameterInformation()
  setPopulationParameterInformation(ka_pop = list(initialValue = 1.31, method = "MLE"),
                                    TVV_pop = list(initialValue = 28.57, method = "MLE"),
                                    TVCl_pop = list(initialValue = 3.52, method = "MLE"),
                                    # gamma_TVV = list(initialValue = gamma_V, method = "MLE"),
                                    omega_TVCl = list(initialValue = omega_Cl, method = "MLE"),
                                    omega_TVV = list(initialValue = omega_V, method = "MLE"),
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
    mutate(param = gsub("TV|_pop", "", param)) %>%
    filter(!param %in% c("a", "b"))
    
  # Get Shrinkage 
     # In Monolix, the mode is calculated via the EBEs task, while the mean is calculated via the Conditional distribution task, as the average of all samples drawn from the conditional distribution.
  shrinkage <- tibble(param = gsub("TV", "", getEtaShrinkage()$conditionalMode$parameters),
                      estimate = getEtaShrinkage()$conditionalMode$shrinkage)

  # # Save the results
  # for (parami in params_for_estimation) {
  #   # Thetas
  #   if (str_detect(parami, "_theta$")) results_matrix[j, parami] <- pull(filter(estimates, param == gsub("_theta$", "", parami)), estimate)
  #   # RSE thetas
  #   if (str_detect(parami, "_theta_rse")) results_matrix[j, parami] <- pull(filter(estimates, param == gsub("_theta_rse", "", parami)), rses)
  #   # Etas
  #   if (str_detect(parami, "_eta$")) results_matrix[j, parami] <- pull(filter(estimates, param == gsub("_eta$", "", parami)), estimate)
  #   # RSE etas
  #   if (str_detect(parami, "_eta_rse")) results_matrix[j, parami] <- pull(filter(estimates, param == gsub("_eta_rse", "", parami)), rses)
  #   # Eta Shrinkage
  #   if (str_detect(parami, "_eta_shr")) results_matrix[j, parami] <-  pull(filter(shrinkage, param == gsub("_eta_shr", "", parami)), estimate)
  # }
  
  # Save the results
  for (parami in params_for_estimation) {
    # Determine the action based on the suffix of `parami`
    action <- case_when(
      str_detect(parami, "_theta$") ~ "theta",
      str_detect(parami, "_theta_rse") ~ "theta_rse",
      str_detect(parami, "_eta$") ~ "eta",
      str_detect(parami, "_eta_rse") ~ "eta_rse",
      str_detect(parami, "_eta_shr") ~ "eta_shr",
      TRUE ~ NA_character_ # Default case
    )
    
    if (!is.na(action)) {
      results_matrix[j, parami] <- switch(
        action,
        "theta" = pull(filter(estimates, param == gsub("_theta$", "", parami)), estimate),
        "theta_rse" = pull(filter(estimates, param == gsub("_theta_rse", "", parami)), rses),
        "eta" = pull(filter(estimates, param == paste0("omega_", gsub("_eta$", "", parami))), estimate),
        "eta_rse" = pull(filter(estimates, param == paste0("omega_", gsub("_eta_rse", "", parami))), rses),
        "eta_shr" = pull(filter(shrinkage, param == gsub("_eta_shr", "", parami)), estimate)
      )
    }
  }
  
  message(paste0("##################\n##################\n##################\n", 
                 round(j * 100/n_datasets, 2), 
                 "% complete\n##################\n##################\n##################"))
}
write.csv(results_matrix, paste0("base_model_refit/refit.csv"), row.names = FALSE, quote = FALSE)

# Obtain plots
# ------------
# Reshape the data for plotting
dat <- read_csv(paste0("base_model_refit/refit.csv"), show_col_types = FALSE)
colnames(dat) <- gsub("_shr", "shr", gsub("_rse", "rse", colnames(dat)))
  
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
         metric = gsub("ETAsSHR", "ETA SHRINKAGE", metric),
         metric = factor(metric, levels = unique(metric)),
         parameter = factor(parameter, levels = unique(parameter))) %>% 
  filter(value < 1000) # Filter out one V ETA RSE that is 1583


# Create the boxplots with facets for metrics (theta, eta, shrinkage)
ggplot(dat_long, aes(x = parameter, y = value, colour = parameter)) +
  geom_boxplot() +
  facet_wrap(~metric, scales = "free", nrow = 3) +
  theme_bw() +
  labs(
    title = "Boxplots of Estimates for Parameters Ka, Cl, and V",
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
ggsave("base_model_estimates.png", height = 12 * x, width = 9 * x)









