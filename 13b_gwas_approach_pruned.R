# Load libraries
# --------------
library(tidyverse)
library(data.table)

sim_date <-"23_12_2024"
n_datasets <- 100
scenarios <- c("low", "high")
n_snps <- c("10_3", "10_4", "10_5", "10_6")
n_true_snps <- 9 
topx <- 100

for (effect_scenario in scenarios) {
  # Add folders if they don't exist
  if (!dir.exists(effect_scenario)) dir.create(effect_scenario)
  folder_path <- paste0(effect_scenario, "/gwas")
  if (!dir.exists(folder_path)) dir.create(folder_path)
  
  for (n_snp in n_snps) {
    
    method_list <- tibble(Cl = vector("character", n_datasets))
    time_list <- tibble(end = vector("character", n_datasets),
                        start = vector("character", n_datasets))
    
    for (j in 1:n_datasets) {
      # Get time of code execution start
      start <- Sys.time()
      
      # Get DV data to get dict
      dict <- read_csv(paste0("sim_DV_", effect_scenario, "/sim_data", j, "_", sim_date, ".csv"), show_col_types = FALSE) %>%
        select(id = ID, ID2) %>%#
        distinct()
      
      # Get ETAs (to use as outcome)
      etas <- read_csv(paste0(effect_scenario, "/base_model/sim_", j, "_", "etas.csv"), show_col_types = FALSE) %>%
        filter(OCC == 1) %>%
        left_join(dict) %>%
        select(FID = ID2, IID = ID2, eta_TVCl) 
      
      # Save etas these as phenotype file
      write.table(etas, paste0("pheno_", n_snp, ".txt"), sep = " ", row.names = FALSE, quote = FALSE, na = "NA")
      
      # Run GWAS in plink
      dat_path <- paste0("pruned_", n_snp, "/dat_", n_snp, "_", j)
      system(paste0("/pub59/iasiimwe/plink2 --bfile ", dat_path, " --pheno ", 
                    paste0("pheno_", n_snp, ".txt --glm dominant --out "),  paste0(effect_scenario, "/gwas/gwas_pruned_",n_snp, "_", j)))

      # Read the results
      results <- as_tibble(fread(paste0(effect_scenario, "/gwas/gwas_pruned_",n_snp, "_", j, ".eta_TVCl.glm.linear"))) 
      colnames(results) <- c("CHR", "BP", "SNP", "alleleA", "alleleB", "A1", "TEST", "OBS_CT", "BETA", "SE", "T_STAT", "P")
      
      # Get true SNPs
      true_snps <- fread(paste0("true_covar/true_snps_", j, ".txt"), header = FALSE, sep = " ") %>%
        rename(SNP = V1) %>%
        mutate(SNP2 = paste0("SNP", row_number()))
      
      # Add to results
      results <- results %>%
        left_join(true_snps) %>%
        mutate(SNP2 = if_else(is.na(SNP2), SNP, SNP2)) %>%
        select(-SNP) %>%
        rename(SNP = SNP2) %>%
        relocate(SNP, .after = BP) %>%
        arrange(P)
      
      end <- Sys.time()
      
      topx_to_use <- min(topx, nrow(results))
      method_list$Cl[j] <- paste0(results$SNP[1:topx_to_use], collapse = "; ")
      time_list$start[j] <- paste("T", as.character(ymd_hms(start)))
      time_list$end[j] <- paste("T", as.character(ymd_hms(end)))
      
      # Save Bonferroni-corrected significant results
      write.csv(results[results$P < 0.05/nrow(results),],  
                file = paste0(effect_scenario, "/gwas/gwas_pruned_",n_snp, "_", j, ".csv"), 
                row.names = FALSE)
      
      # The Manhattan plots can be generated separately, so let us get a plot folder to save all the results
      if (!dir.exists(paste0(effect_scenario, "/gwas/plot"))) dir.create(paste0(effect_scenario, "/gwas/plot"))
      write.csv(results,  
                file = paste0(effect_scenario, "/gwas/plot/gwas_pruned_",n_snp, "_", j, ".csv"), 
                row.names = FALSE)

      message(paste0(round(j*100/n_datasets, 3), "% complete"))
    }
    write.csv(time_list, paste0(folder_path, "/pruned_", n_snp, "_time_gwas.csv"), row.names = FALSE)
    write.csv(method_list, paste0(folder_path, "/pruned_", n_snp, "_results_gwas.csv"), row.names = FALSE) 
  }
}
