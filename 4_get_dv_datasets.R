# Load required packages
# ---------------------
library(lixoftConnectors)
initializeLixoftConnectors(software = "simulx", path = "/pub59/iasiimwe/Lixoft/MonolixSuite2024R1/", force = TRUE)
library(tidyverse)
library(data.table)
library(rxode2)
library(MASS)
library(ggh4x)

select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate

# Add folders if they don't exist
if (!dir.exists("sim_DV")) dir.create("sim_DV")

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

# Without re-inflation
omega_ka <- omega_squared_fn(68.5)
omega_Cl <- omega_squared_fn(17.5)
omega_V <- omega_squared_fn(18.2)
gamma_V <- omega_squared_fn(11.7)

sim_date <-"23_12_2024"

# Load clinical covariates dataset
clin_dat <- read_csv(paste0("clin_cov_dataset_", sim_date, ".csv"), show_col_types = FALSE) %>%
  rename(ID2 = ID) %>%
  mutate(ID = row_number(), .after = ID2)

# Time_tb
time_per_day <- c(0, 0.3, 0.9, 2.2, 4.5, 8)
time_tb <- tibble(TIME = rep(time_per_day, 2),
                  EVID = rep(c(1, c(rep(0, length(time_per_day) - 1))), 2),
                  OCC = c(rep(1, length(time_per_day)), rep(2, length(time_per_day))),
                  DV = c(".", 8.46, 24.44, 56.39, 38.81, 30.81, ".", 4.27, 33.27, 40.66, 32.77, 25.31)) # From SUBID 1500
dv_tb <- tibble()
for (i in 1:nrow(clin_dat)) dv_tb <- dv_tb %>% bind_rows(mutate(time_tb, ID = i))

set.seed(7) # For reproducibility
time_sd <- 0
effect_size <- 0.5

n_datasets <- 100
for (j in 1:n_datasets) {
  # Convert to ped format
  system(paste0("C:/Users/iasiimww/Downloads/plink/plink --bfile datasets/10_3/dat_10_3_", j, " --recode tab --out temp"))
  
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
  
  # Join to the clinical data
  clin_snp_dat <- clin_dat %>%
    left_join(ped_snps) %>%
    relocate(ID2, .after = SNP9)
  
  # Add 'DV' data
  dat <- clin_snp_dat %>%
    left_join(dv_tb) %>%
    relocate(TIME:DV, .after = ID) %>%
    mutate(AMT = round(WT * 25/400) * 400, # 25 mg/kg, rounded to the nearest 400 mg
           AMT = ifelse(TIME == 0 & EVID == 1, AMT, 0), .after = OCC) 
  write.csv(select(dat, -ID2), "temp.csv", row.names = FALSE)
  
  # Start simulation
  initializeLixoftConnectors(software = "monolix", path = "/pub59/iasiimwe/Lixoft/MonolixSuite2024R1/", force = TRUE)
  newProject(data = list(dataFile = 'temp.csv',
                         headerTypes = c("id", "time", "evid", "occ", "amount", "observation", 
                                         "catcov", "regressor", rep("catcov", 9))),
             modelFile = 'tb_base_vinnard.txt')
  
  # Change error model, distribution, variability and correlations
  setErrorModel(DV = "combined1")
  setIndividualParameterDistribution(ka = "logNormal", TVV = "logNormal", TVCl = "logNormal")
  setIndividualParameterVariability(id = list (ka = TRUE, TVV = TRUE, TVCl = TRUE),
                                    OCC = list (ka = FALSE, TVV = TRUE, TVCl = FALSE))
  # getIndividualParameterModel() to confirm changes
  
  # Add sex and SNPs as covariates
  setCovariateModel(TVCl = c(
    #SEX = TRUE, 
                             SNP1 = TRUE, SNP2 = TRUE, SNP3 = TRUE, SNP4 = TRUE, 
                             SNP5 = TRUE, SNP6 = TRUE, SNP7 = TRUE, SNP8 = TRUE, SNP9 = TRUE))
  
  # getPopulationParameterInformation()
  setPopulationParameterInformation(ka_pop = list(initialValue = 1.31, method = "FIXED"),
                                    TVV_pop = list(initialValue = 28.57, method = "FIXED"),
                                    TVCl_pop = list(initialValue = 3.52, method = "FIXED"),
                                    #beta_TVCl_SEX_2 = list(initialValue = -0.4, method = "FIXED"),
                                    beta_TVCl_SNP1_1 = list(initialValue = effect_size, method = "FIXED"),
                                    beta_TVCl_SNP2_1 = list(initialValue = effect_size, method = "FIXED"),
                                    beta_TVCl_SNP3_1 = list(initialValue = effect_size, method = "FIXED"),
                                    beta_TVCl_SNP4_1 = list(initialValue = effect_size, method = "FIXED"),
                                    beta_TVCl_SNP5_1 = list(initialValue = effect_size, method = "FIXED"),
                                    beta_TVCl_SNP6_1 = list(initialValue = effect_size, method = "FIXED"),
                                    beta_TVCl_SNP7_1 = list(initialValue = effect_size, method = "FIXED"),
                                    beta_TVCl_SNP8_1 = list(initialValue = effect_size, method = "FIXED"),
                                    beta_TVCl_SNP9_1 = list(initialValue = effect_size, method = "FIXED"),
                                    gamma_TVV = list(initialValue = gamma_V, method = "FIXED"),
                                    omega_TVCl = list(initialValue = omega_Cl, method = "FIXED"),
                                    omega_TVV = list(initialValue = omega_V, method = "FIXED"),
                                    omega_ka = list(initialValue = omega_ka, method = "FIXED"),
                                    a = list(initialValue = 2.41, method = "FIXED"),
                                    b = list(initialValue = 0.22, method = "FIXED"),
                                    c = list(initialValue = 1, method = "FIXED"))
  
  # Set tasks in scenario
  scenario <- getScenario()
  scenario$tasks = c(populationParameterEstimation = TRUE, 
                     conditionalDistributionSampling = TRUE,
                     conditionalModeEstimation = TRUE, 
                     standardErrorEstimation = FALSE, 
                     logLikelihoodEstimation = FALSE,
                     plots = FALSE)
  scenario$linearization = FALSE
  setScenario(scenario)
  
  runScenario()
  # getPopulationParameterInformation()
  saveProject("temp.mlxtran")
  
  initializeLixoftConnectors(software = "simulx", path = "/pub59/iasiimwe/Lixoft/MonolixSuite2024R1/", force = TRUE)
  importProject("temp.mlxtran")
  setPreferences(exportsimulationfiles = TRUE)
  runSimulation()
  sim_results <- getSimulationResults()$res$DV %>%
    select(ID = original_id, OCC, TIME = time, DV) %>%
    mutate(ID = as.numeric(ID))
  sim_results <- select(dat, -DV) %>%
    left_join(sim_results) %>%
    mutate(NTIME = TIME, .after = ID) %>%
    relocate(ID2, .after = DV) %>%
    mutate(DV = ifelse(is.na(DV), ".", as.character(round(DV, 3)))) %>%
    relocate(DV, .after = AMT) %>%
    mutate(TIME = round(rnorm(nrow(.), NTIME, time_sd), 3), # Add some variability if necessary
           TIME = ifelse(NTIME == 0, 0, TIME)) 
  
  if (TRUE) {
    subject_ids <- unique(sim_results$ID)
    to_include <- sample(subject_ids, ceiling(0.6 * length(subject_ids))) # We only want 60% subjects for occasion 2
    to_exclude <- subject_ids[!subject_ids %in% to_include]
    sim_results <- sim_results %>% filter(!(ID %in% to_exclude & OCC == 2))
  }
  
  write.csv(sim_results, paste0("sim_DV/sim_data", j, "_", sim_date, ".csv"), row.names = FALSE, quote = FALSE)
  
  message(paste0("##################\n##################\n##################\n", 
                 round(j * 100/n_datasets, 2), 
                 "% complete\n##################\n##################\n##################"))
}

# Check the first four datasets
nonmem_tb <- tibble()
for (j in 1:4) nonmem_tb <- bind_rows(nonmem_tb, mutate(read_csv(paste0("sim_DV/sim_data", j, "_", sim_date, ".csv"), show_col_types = FALSE), Simulation = paste0("Dataset ", j)))

# Create the concentration-time plots
nonmem_tb %>%
  mutate(DV = as.numeric(if_else(DV == ".", NA, DV)),
         OCC = ifelse(OCC == 1, "After 7 treatment days", "After 56 treatment days"),
         OCC = factor(OCC, levels = unique(OCC))) %>%
  ggplot(aes(TIME, DV, color = OCC, group = interaction(ID, OCC))) +
  geom_line() +    # Connect points with lines
  geom_point() +   # Add points to the plot
  facet_nested_wrap(vars(Simulation, OCC), nrow = 2) + # , scales = "free_y"
  labs(
    # title = "Time-Concentration Plot by Subject and Occasion",
    x = "Time after dose (hours)",
    y = "Concentration (mg/L)",
    color = "Occasion"
  ) +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0, vjust = 0, size = 14, face = "bold"))
x <- 0.8
ggsave("simulated_dv_snps.png", height = 9 * x, width = 16 * x)

