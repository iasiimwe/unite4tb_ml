# Load libraries
# --------------
library(lixoftConnectors)
# initializeLixoftConnectors(software = "monolix", path = 'C:/Program Files/Lixoft/MonolixSuite2024R1', force = TRUE)
initializeLixoftConnectors(software = "monolix", path = "/pub59/iasiimwe/Lixoft/MonolixSuite2023R1/", force = TRUE)
library(tidyverse)
library(data.table)

# Slightly modify the model building function
runModelBuilding2 <- function (...) {
  arguments = list(...)
  name = names(arguments)
  offset = grep("settings", name)
  if (length(offset) != 0) {
    inputs = arguments[[offset]]
    settings = list()
    inputNames = names(inputs)
    for (i in 1:length(inputs)) {
      settings[[inputNames[i]]] = inputs[[i]]
    }
    if (!is.null(settings$useLin) && !is.logical(settings$useLin)) {
      .error("Unexpected type encountered for \"useLin\" field. Please give a boolean.")
      return(invisible(FALSE))
    }
    if (length(settings$parameters) > 0) 
      settings$parameters = as.list(settings$parameters)
    if (length(settings$covariates) > 0) 
      settings$covariates = as.list(settings$covariates)
    if (!is.null(settings$threshold) && !is.null(settings$threshold$correlation)) {
      if (!is.numeric(settings$threshold$correlation) || 
          length(settings$threshold$correlation) != 2) {
        .error("Unexpected type encountered for \"correlation\" field. Please give 2 numeric values")
        return(invisible(FALSE))
      }
      if (min(settings$threshold$correlation <= 0) || max(settings$threshold$correlation) >= 
          1) {
        .error("Unexpected value(s) encountered for \"correlation\" field. Please give 2 numeric values stricly between 0 and 1.")
        return(invisible(FALSE))
      }
    }
    if (!is.null(settings$threshold) && !is.null(settings$threshold$lrt)) {
      if (!is.numeric(settings$threshold$lrt) || length(settings$threshold$lrt) != 
          2) {
        .error("Unexpected type encountered for \"lrt\" field. Please give 2 numeric values")
        return(invisible(FALSE))
      }
      if (min(settings$threshold$lrt <= 0) || max(settings$threshold$lrt) >= 
          1) {
        .error("Unexpected value(s) encountered for \"lrt\" field. Please give 2 numeric values stricly between 0 and 1.")
        return(invisible(FALSE))
      }
    }
  }
  else {
    settings = getModelBuildingSettings()
  }
  output = .processRequest(.software(), "setmodelbuildingarguments", 
                           settings, "synchronous", type = "STATUS")
  if (output) 
    output = .processRequest(.software(), "runmodelbuilding", 
                             list(), "synchronous", type = "STATUS")
  return(invisible(output))
}
environment(runModelBuilding2) <- environment(runModelBuilding)

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

## Random effects (considering log distribution)
omega_squared_fn <- function(percent_CV) log((percent_CV/100)^2 + 1)

# # Without re-inflation
omega_ka <- omega_squared_fn(68.5)
omega_Cl <- omega_squared_fn(17.5)

sim_date <-"23_12_2024"
effect_scenario <- c("high")
n_snp <- c("10_3")
j <- 1 # dataset 1
p_forwards <- 0.05 # forward p-value threshold for LRT
p_backwards <- 0.01

ml_methods <- c("scm", "cossac", "cossac") # The third will be "covSAMBA-COSSAC"
ml_method_names <- c("scm", "cossac", "cossac_samba")
use_cossac <- c(FALSE, FALSE, TRUE)
k <- 3

# Get DV data
dv_dat <- read_csv(paste0("sim_DV_", effect_scenario, "/sim_data", j, "_", sim_date, ".csv"), show_col_types = FALSE) %>%
  select(!c(contains("SNP"), NTIME))

# Get SNP data (convert to ped format, and process)
dat_path <- paste0(n_snp, "/dat_",n_snp, "_", j)
temp_path <- paste0("cossac_temp_",n_snp)
# system(paste0("C:/Users/asiimwe/Downloads/plink_win64_20241022/plink.exe --bfile ", dat_path, " --recode tab --out ", temp_path), intern = TRUE)
system(paste0("/pub59/iasiimwe/plink1.9/plink --bfile ", dat_path, " --recode tab --out ", temp_path))
ped <- tibble(fread(paste0(temp_path, ".ped"), header = FALSE))
map <- tibble(fread(paste0(temp_path, ".map"), header = FALSE))
snps <- pull(map, V2)
ped_start <- ped %>% select(V1:V6) # This selects FID, IID, father's ID, mother's ID, sex and phenotype data
snp_columns <- colnames(ped)[!colnames(ped) %in% colnames(ped_start)]
ped_snps <- ped %>%
  select(all_of(snp_columns)) 
colnames(ped_snps) <- snps 

# True covariates - we need to put them as the first 9 covariates for later tracking
true_snps <- fread(paste0("true_covar/true_snps_", j, ".txt"), header = FALSE, sep = " ")$V1
ped_snps <- ped_snps %>%
  select(any_of(true_snps), everything()) %>% # This allows the order to be the same as in 'true_snps'
  mutate_all(dominant_coding_fn) %>% # change to 'dominant' coding
  mutate(ID2 = ped_start$V1)
colnames(ped_snps) <- c(paste0("SNP", 1:(ncol(ped_snps) - 1)), "ID2")
ped_snps <- ped_snps %>%
  relocate(ID2, .before = 1)

# Let us get the first 100 SNPs and join them to the dv data
n_snps <- seq(0, 500, by = 50)

dat_all <- ped_snps[, 1:(max(n_snps) + 1)] %>%
  left_join(dv_dat) %>%
  relocate(ID:WT, .before = 1) %>%
  select(-ID2)

time_tb <- tibble(n_snps = n_snps,
                  end = vector("character", length(n_snps)),
                  start = vector("character", length(n_snps)))

for (i in seq_along(n_snps)) {
  message(paste0("Starting ", n_snps[[i]], " SNPs"))
  
  # Get time of code execution start
  start <- Sys.time()
  
  ncols <- 8 + n_snps[[i]]
  
  dat <- dat_all[, 1:ncols]
  write.csv(dat, "temp_cossac.csv", row.names = FALSE)
  
  header_types <- c("id", "time", "evid", "occ", "amount", "observation", "catcov", "regressor",
                    rep("catcov", n_snps[[i]])) 
  
  initializeLixoftConnectors(software = "monolix", path = "/pub59/iasiimwe/Lixoft/MonolixSuite2023R1/", force = TRUE)
  newProject(data = list(dataFile = "temp_cossac.csv",
                         headerTypes = header_types),
             modelFile = 'tb_base_vinnard.txt')
  
  # Change error model, distribution, variability and correlations
  setErrorModel(DV = "combined1")
  setIndividualParameterDistribution(ka = "logNormal", TVV = "logNormal", TVCl = "logNormal")
  setIndividualParameterVariability(id = list (ka = TRUE, TVV = FALSE, TVCl = TRUE),
                                    OCC = list (ka = FALSE, TVV = FALSE, TVCl = FALSE))
  # getIndividualParameterModel() to confirm changes
  
  setPopulationParameterInformation(ka_pop = list(initialValue = 1.31, method = "MLE"),
                                    TVV_pop = list(initialValue = 28.57, method = "FIXED"),
                                    TVCl_pop = list(initialValue = 3.52, method = "MLE"),
                                    omega_TVCl = list(initialValue = omega_Cl, method = "MLE"),
                                    omega_ka = list(initialValue = omega_ka, method = "MLE"),
                                    a = list(initialValue = 2.41, method = "FIXED"),
                                    b = list(initialValue = 0.22, method = "MLE"),
                                    c = list(initialValue = 1, method = "FIXED"))
  
  # Set tasks in scenario
  scenario <- getScenario()
  scenario$tasks = c(populationParameterEstimation = TRUE, 
                     conditionalDistributionSampling = TRUE,
                     conditionalModeEstimation = FALSE, 
                     standardErrorEstimation = FALSE, 
                     logLikelihoodEstimation = FALSE,
                     plots = FALSE)
  scenario$linearization = FALSE
  setScenario(scenario)
  
  # get the default covariate search settings
  set <- getModelBuildingSettings()
  
  # modify the settings
  set$strategy <- ml_methods[[k]]
  set$useSambaBeforeCossac <- use_cossac[[k]]
  set$criterion <- "LRT" 
  set$threshold$lrt[1] <- p_forwards
  set$threshold$lrt[2] <- p_backwards
  set$parameters <- c("TVCl") # Only on Cl
  
  # No linearization
  set$useLin <- FALSE
  
  # run the covariate search
  runModelBuilding2(settings = set)
  
  # get the results. The best model is indicated with $bestModel = T
  mod_results <- getModelBuildingResults()
  # str(mod_results)
  
  
  best_mod <- c()
  for (m in seq_along(mod_results)) if (mod_results[[m]]$bestModel == TRUE) best_mod <- c(best_mod, m)
  
  # method_list[[k]]$Cl[j] <- mod_results[[best_mod]]$individualModels["TVCl", ] %>%
  #   tibble() %>%
  #   mutate(x = NA) %>%
  #   pivot_longer(-x) %>%
  #   filter(value == TRUE) %>%
  #   pull(name) %>%
  #   paste0(., collapse = ", ")
  
  end <- Sys.time()
  time_tb$start[i] <- paste("T", as.character(ymd_hms(start)))
  time_tb$end[i] <- paste("T", as.character(ymd_hms(end)))
}
write_csv(time_tb, "cossac_test.csv")


time_tb <- read_csv("cossac_test.csv", show_col_types = FALSE) %>%
  mutate(end = ymd_hms(gsub("T ", "", end)),
         start = ymd_hms(gsub("T ", "", start)),
         time_h = as.numeric(end-start, "hours")) %>%
  select(`Number of SNPs` = n_snps, `Time (hours)` = time_h)


ggplot(time_tb, aes(x = `Number of SNPs`, y = `Time (hours)`)) +
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ poly(x, 3), se = TRUE, color = "blue") +
  ggtitle("Modeling Runtime of COSSAC-SAMBA by Number of SNPs") +
  theme_bw() + 
  theme(axis.text = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0, vjust = 0, size = 12, face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  # annotate("text", x = 400, y = 11, label = "Time == (NSNPs)^3", parse = TRUE, colour = "blue")
  annotate("text", x = 400, y = 11, label = "y == x^3", parse = TRUE, colour = "blue")
x <- 0.6
ggsave("COSSAC_SAMBA_Time.png", height = 10 * x, width = 10 * x)
