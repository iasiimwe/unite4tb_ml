# Load libraries
library(tidyverse)
library(data.table)
library(RColorBrewer)
library(patchwork)
library(export)
library(qqman) 

# Relevant functions
dominant_coding_fn <- function(x) {
  x <- gsub(" ", "", x) 
  x[x == "00"] <- NA
  y <- table(x)
  y <- names(sort(-y)) # Sort in descending order, assume mutant type is the least abundant
  if (length(y) == 3) x <- gsub(y[[3]], "1", x) # Give mutant homozygotes "2"
  if (length(y) > 1) x <- gsub(y[[2]], "1", x) # Give heterozygotes "1" 
  x <- gsub(y[[1]], "0", x) # Give homozygotes 0
  return(as.numeric(x))
}

# Get datasets
n_datasets <- 100

mechanisms <- c("MCAR", "MAR", "MNAR")
missing_percentages <- c(5, 10, 20, 50)
set.seed(7)

results_tb <- tibble()
for (mechanism in mechanisms) {
  for (missing_percentage in missing_percentages) {
    message(paste0("Starting\n     Mechanism: ", 
                   mechanism, "\n     missing %: ", missing_percentage))
    
    for (j in 1:n_datasets) {
      # True SNPs
      true_snps <- fread(paste0("true_covar/true_snps_", j, ".txt"), header = FALSE, sep = " ")$V1
      
      # Get original dataset
      dat_path <- paste0("10_3/dat_10_3_", j)
      system(paste0("/pub59/iasiimwe/plink1.9/plink --bfile ", dat_path, " --recode tab --out temp_fig1"))
      ped <- tibble(fread("temp_fig1.ped", header = FALSE))
      map <- tibble(fread("temp_fig1.map", header = FALSE))
      snps <- pull(map, V2)
      ped_start <- ped %>% select(V1:V6) # This selects FID, IID, father's ID, mother's ID, sex and phenotype data
      snp_columns <- colnames(ped)[!colnames(ped) %in% colnames(ped_start)]
      ped_snps <- ped %>%
        select(all_of(snp_columns)) 
      colnames(ped_snps) <- snps 
      ped_snps <- ped_snps %>%
        select(any_of(true_snps)) %>% # This allows the order to be the same as in 'true_snps'
        mutate_all(dominant_coding_fn) %>%
        mutate(ID2 = ped_start$V1, Sex = ped_start$V5)
      colnames(ped_snps) <- c(paste0("OSNP", 1:(ncol(ped_snps) - 2)), "ID2", "Sex")
      ped_snpsO <- ped_snps
      
      # Get missing dataset
      dat_path <- paste0("missing/", missing_percentage, "/dat_", mechanism, "_", j)
      system(paste0("/pub59/iasiimwe/plink1.9/plink --bfile ", dat_path, " --recode tab --out temp_fig1"))
      ped <- tibble(fread("temp_fig1.ped", header = FALSE))
      map <- tibble(fread("temp_fig1.map", header = FALSE))
      snps <- pull(map, V2)
      ped_start <- ped %>% select(V1:V6) # This selects FID, IID, father's ID, mother's ID, sex and phenotype data
      snp_columns <- colnames(ped)[!colnames(ped) %in% colnames(ped_start)]
      ped_snps <- ped %>%
        select(all_of(snp_columns)) 
      colnames(ped_snps) <- snps 
      ped_snps <- ped_snps %>%
        select(any_of(true_snps)) %>% # This allows the order to be the same as in 'true_snps'
        mutate_all(dominant_coding_fn) %>%
        mutate(ID2 = ped_start$V1, Sex = ped_start$V5)
      colnames(ped_snps) <- c(paste0("MSNP", 1:(ncol(ped_snps) - 2)), "ID2", "Sex")    

      # Combine the datasets
      dat <- ped_snpsO %>%
        left_join(ped_snps) %>%
        mutate(Mechanism = mechanism,
               Missing_percentage = missing_percentage,
               Dataset = j)
      results_tb <- bind_rows(results_tb, dat)
      
      message(paste0("##################\n##################\n##################\n", 
                     round(j * 100/n_datasets, 2), 
                     "% complete\n     Mechanism: ", 
                     mechanism, "\n     missing %: ", missing_percentage, 
                     "\n##################\n##################\n##################"))
    }
  }
}

dat_all <- tibble()
for (mechanism in unique(results_tb$Mechanism)) {
  for (missing_percentage in unique(results_tb$Missing_percentage)) {
    for (dataset in unique(results_tb$Dataset)) {
      message(paste0("Starting\n     Mechanism: ", mechanism, "\n     Missing %: ", 
                     missing_percentage, "\n     Dataset: ", dataset))
      dat <- results_tb %>%
        filter(Mechanism == mechanism, Missing_percentage == missing_percentage, Dataset == dataset)
      # Get long data
      ID_df <- distinct(dat, ID2, Sex, Mechanism, Missing_percentage, Dataset) %>%
        mutate(Sex = factor(if_else(Sex == 1, "Male", "Female")))
      original_SNPs_df <- dat %>%
        select(ID2, starts_with("OSNP")) %>%
        pivot_longer(cols = starts_with("OSNP"), names_to = "SNP", values_to = "Original") %>%
        mutate(SNP = gsub("O", "", SNP))
      missing_SNPs_df <- dat %>%
        select(ID2, starts_with("MSNP")) %>%
        pivot_longer(cols = starts_with("MSNP"), names_to = "SNP", values_to = "Missing") %>%
        mutate(SNP = gsub("M", "", SNP))
      
      dat_long <- tibble(ID_df) %>%
        left_join(original_SNPs_df) %>%
        left_join(missing_SNPs_df) %>%
        rename(ID = ID2) %>%
        mutate(Original = factor(if_else(Original == 0, "Wild-type", "Mutant"), levels = c("Wild-type", "Mutant")),
               Missing = if_else(is.na(Missing), NA, 1))
      
      dat_all <- dat_all %>%
        bind_rows(dat_long)
    }
  }
}
dat_all

write_csv(dat_all, "fig_1.csv")
system("/pub59/iasiimwe/htslib/htslib-1.3.2/bgzip -f fig_1.csv")


# Missingness mechanisms
results_tb <- tibble(fread("fig_1.csv.gz"))
dat <- results_tb %>%
  group_by(Sex, Mechanism, Original) %>%
  summarise(Proportion_NA = mean(is.na(Missing)) * 100, .groups = "drop") %>%
  mutate(Sex = factor(Sex),
         Mechanism = factor(Mechanism, levels = c("MCAR", "MAR", "MNAR")),
         # Allele = recode(Original, "Mutant" = "≥1 Mutant"),
         `SNP-type` = factor(Original, levels = c("Wild-type", "Mutant")))
pl <- ggplot(dat, aes(x = Mechanism, y = Proportion_NA, fill =`SNP-type`)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Sex) +
  labs(x = "Mechanism", y = "Percentage of NAs",
       title = "Missingness by Mechanism and SNP-type") +
  scale_fill_manual(values = c("Wild-type" = "#7CAE00", "Mutant" =  "#C77CFF")) + # c("Wild-type" = "#00BFC4", "Mutant" = "#F8766D")
  theme_bw() +
  theme(legend.position = "top",
        strip.text = element_text(size = 12),
        axis.text = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 12),
        # axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0),
        legend.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5, vjust = 0, size = 12, face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
x <- 0.3
graph2ppt(pl, "Fig_1", append = FALSE, width = 16 * x, height = 9 * x)


# Vinnard original DV data
pl <- read_csv("S1Table.csv", show_col_types = FALSE) %>%
  distinct(ID = SUBID, OCC, TIME, DV = CONC) %>%
  mutate(OCC = ifelse(OCC == 1, "After ~20 treatment days", "After ~74 treatment days"),
         OCC = factor(OCC, levels = unique(OCC))) %>%
  ggplot(aes(TIME, DV, color = OCC, group = interaction(ID, OCC))) +
  facet_wrap(~OCC, nrow = 1) + 
  geom_line() +    
  geom_point() +   
  labs(title = "Vinnard et al. Pyrazinamide PK dataset", x = "Time after dose (hours)", y = "Concentration (mg/L)", color = "Occasion") +
  theme_bw() +
  theme(legend.position = "none",
        strip.text = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, vjust = 0, size = 12, face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
graph2ppt(pl, "Fig_1", append = TRUE, width = 16 * x, height = 9 * x)

# Vinnard covariate dataset
dat <- read_csv("S1Table.csv", show_col_types = FALSE) %>%
  distinct(SUBID, SEX = GENDER, WT) %>% # based on the WT and HT means, it is likely that 0 == F, and 1 == M
  group_by(SUBID) %>%
  slice_head(n = 1) 
# Number of males
n_male <- dat %>% filter(SEX == 1) %>% nrow()
n_female <- dat %>% filter(SEX == 2) %>% nrow()
# Calculate medians and correlation
median_wt <- median(dat$WT)
# Calculate sex-specific medians
sex_medians <- dat %>%
  group_by(SEX) %>%
  summarise(WT = median(WT), .groups = "drop")
male_wt <- sex_medians %>% filter(SEX == 1) %>% pull(WT)
female_wt <- sex_medians %>% filter(SEX == 2) %>% pull(WT)
# Weight ratio
ratio_WT <- male_wt/female_wt
# Create density plot for sex-specific weight distribution
pl <- ggplot(dat, aes(x = WT, color = factor(SEX), fill = factor(SEX))) +
  geom_density(alpha = 0.4) +
  scale_color_manual(values = c("blue", "pink"), 
                     labels = c(paste0("Male (N = ", n_male, ")"), paste0("Female (N = ", n_female, ")")), 
                     name = "Sex") +
  scale_fill_manual(values = c("blue", "pink"), 
                    labels = c(paste0("Male (N = ", n_male, ")"), paste0("Female (N = ", n_female, ")")), 
                    name = "Sex") +
  labs(title = "Sex-Specific Weight Density Distribution",
       x = "Weight (WT, kg)", y = "Density", caption = sprintf("Overall Median WT: %.0f", median_wt)) +
  theme_bw() +
  theme(plot.caption = element_text(hjust = 0.5, size = 10),
        legend.position = "top",
        legend.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, vjust = 0, size = 12, face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
graph2ppt(pl, "Fig_1", append = TRUE, width = 16 * x, height = 9 * x)


# Simulated covariate dataset
sim_date <-"23_12_2024"
dat <- read_csv(paste0("clin_cov_dataset_", sim_date, ".csv"), show_col_types = FALSE)
n_male <- dat %>% filter(SEX == 1) %>% nrow()
n_female <- dat %>% filter(SEX == 2) %>% nrow()
median_wt <- median(dat$WT)
sex_medians <- dat %>%
  group_by(SEX) %>%
  summarise(WT = median(WT), .groups = "drop")
male_wt <- sex_medians %>% filter(SEX == 1) %>% pull(WT)
female_wt <- sex_medians %>% filter(SEX == 2) %>% pull(WT)
ratio_WT <- male_wt/female_wt
pl <- ggplot(dat, aes(x = WT, color = factor(SEX), fill = factor(SEX))) +
  geom_density(alpha = 0.4) +
  scale_color_manual(values = c("blue", "pink"), 
                     labels = c(paste0("Male (N = ", n_male, ")"), paste0("Female (N = ", n_female, ")")),  
                     name = "Sex") +
  scale_fill_manual(values = c("blue", "pink"), 
                    labels = c(paste0("Male (N = ", n_male, ")"), paste0("Female (N = ", n_female, ")")), 
                    name = "Sex") +
  labs(title = "Sex-Specific Weight Density Distribution", x = "Weight (WT, kg)", y = "Density", caption = sprintf("Overall Median WT: %.0f", median_wt)) +
  theme_bw() +
  theme(plot.caption = element_text(hjust = 0.5, size = 10),
        legend.position = "top",
        legend.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, vjust = 0, size = 12, face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
graph2ppt(pl, "Fig_1", append = TRUE, width = 16 * x, height = 9 * x)


# Legend for missingness
legend_data <- tibble(Percentage = c("5%", "10%", "20%", "50%", "100%"),
                      x = 1:5,  # Positions of the columns
                      y = 1,    # Single row
                      color = c(brewer.pal(7, "Greens")[7], 
                                brewer.pal(7, "Greens")[4], 
                                brewer.pal(7, "Reds")[4], 
                                brewer.pal(7, "Reds")[7],
                                "red"))
pl <- ggplot(legend_data, aes(x = x, y = y, fill = color)) +
  geom_tile() +  # Draw tiles (heatmap blocks)
  geom_text(aes(label = Percentage), vjust = -1, size = 5, colour = "white") +  # Add labels above the tiles
  scale_fill_identity() +  # Use the colors directly without mapping them to a scale
  theme_void() +  # Remove axes and grid lines
  theme(legend.position = "none")  # Remove default legends
graph2ppt(pl, "Fig_1", append = TRUE, width = 16 * x, height = 9 * x)


# Legend for the 9 SNPs
SNPs <- c("SNP1", "SNP2", "SNP3", "SNP4", "SNP5", "SNP6", "SNP7", "SNP8", "SNP9")
MAF <- rep(c("5%", "10%", "20%"), times = 3)  # Repeated for each column
Correlation <- rep(c("Low", "Moderate", "High"), each = 3)  # Repeated for each row
green_colours <- brewer.pal(5, "Greens")
blue_colours <- brewer.pal(5, "Blues")
red_colours <- brewer.pal(5, "Reds")
manual_colours <- c(SNP9 = green_colours[5], SNP8 = green_colours[3], SNP7 = green_colours[2],
                    SNP6 = blue_colours[5], SNP5 = blue_colours[3], SNP4 = blue_colours[2],
                    SNP3 = red_colours[5], SNP2 = red_colours[3], SNP1 = red_colours[2])
df <- tibble(SNPs, MAF, Correlation) %>%
  mutate(MAF = factor(MAF, levels = unique(MAF)),
         Correlation = factor(Correlation, levels = unique(Correlation)))
pl <- ggplot(df, aes(x = MAF, y = Correlation, fill = SNPs)) +
  geom_tile(color = "black") +  # Black grid lines for better visibility
  geom_text(aes(label = SNPs), color = "white", size = 5) +  # Label each cell
  scale_fill_manual(values = manual_colours) +  # Apply manual colors
  scale_x_discrete(position = "top") +  # Move x-axis labels to the top
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        panel.grid = element_blank(),
        legend.position = "none")  # Remove grid lines
graph2ppt(pl, "Fig_1", append = TRUE, width = 16 * x, height = 9 * x)

# Simulated DV data (low effects)
pl <- read_csv("low_sim_data1_23_12_2024.csv", show_col_types = FALSE) %>%
  distinct(ID, OCC, TIME, DV) %>%
  filter(DV != ".") %>%
  mutate(DV = as.numeric(DV),
         OCC = ifelse(OCC == 1, "After ~20 treatment days", "After ~74 treatment days"),
         OCC = factor(OCC, levels = unique(OCC))) %>%
  ggplot(aes(TIME, DV, color = OCC, group = interaction(ID, OCC))) +
  facet_wrap(~OCC, nrow = 1) + 
  geom_line() +    
  geom_point() +   
  labs(title = "Simulated PK dataset 1 (low effect size, 0.15)", x = "Time after dose (hours)", y = "Concentration (mg/L)", color = "Occasion") +
  theme_bw() +
  theme(legend.position = "none",
        strip.text = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, vjust = 0, size = 12, face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
graph2ppt(pl, "Fig_1", append = TRUE, width = 16 * x, height = 9 * x)


# Simulated DV data (high effects)
pl <- read_csv("high_sim_data1_23_12_2024.csv", show_col_types = FALSE) %>%
  distinct(ID, OCC, TIME, DV) %>%
  filter(DV != ".") %>%
  mutate(DV = as.numeric(DV),
         OCC = ifelse(OCC == 1, "After ~20 treatment days", "After ~74 treatment days"),
         OCC = factor(OCC, levels = unique(OCC))) %>%
  ggplot(aes(TIME, DV, color = OCC, group = interaction(ID, OCC))) +
  facet_wrap(~OCC, nrow = 1) + 
  geom_line() +    
  geom_point() +   
  labs(title = "Simulated PK dataset 1 (high effect size, 0.5)", x = "Time after dose (hours)", y = "Concentration (mg/L)", color = "Occasion") +
  theme_bw() +
  theme(legend.position = "none",
        strip.text = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, vjust = 0, size = 12, face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
graph2ppt(pl, "Fig_1", append = TRUE, width = 16 * x, height = 9 * x)



# Model fitting - base model
# --------------------------
dat <- read_csv(paste0("base_model_refit.csv"), show_col_types = FALSE)
colnames(dat) <- gsub("_shr", "shr", gsub("_rse", "rse", colnames(dat)))
dat_long <- dat %>%
  pivot_longer(cols = everything(),
               names_to = c("parameter", "metric"),
               names_sep = "_",
               values_to = "value") %>%
  mutate(metric = toupper(metric),
         metric = gsub("ETA", "ETAs", metric),
         metric = gsub("RSE", " (RSEs)", metric),
         metric = gsub("ETAsSHR", "SHRINKAGE", metric),
         metric = gsub("ka", "Ka", metric),
         metric = factor(metric, levels = unique(metric)),
         parameter = gsub("V", "Vd", parameter),
         parameter = factor(parameter, levels = unique(parameter))) %>% 
  filter(value > 0) %>% # Filter out one outlier
  filter(!str_detect(metric, "RSEs"))
pl <- ggplot(dat_long, aes(x = parameter, y = value, colour = parameter)) +
  geom_boxplot() +
  facet_wrap(~metric, scales = "free", nrow = 1) +
  theme_bw() +
  labs(title = "Base Model (Allometric Scaling)", x = "Parameter", y = "Value") +
  theme(legend.position = "none",
        strip.text = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, vjust = 0, size = 12, face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
graph2ppt(pl, "Fig_1", append = TRUE, width = 18 * x, height = 9 * x)


# Obtain plots
# ------------
# Reshape the data for plotting
dat <- read_csv("with_dv.csv", show_col_types = FALSE)
colnames(dat) <- gsub("b", "b_b", gsub("_rse", "rse", colnames(dat)))
dat_long <- dat %>%
  pivot_longer(cols = everything(),
               names_to = c("parameter", "metric"),
               names_sep = "_",
               values_to = "value") %>%
  mutate(metric = toupper(metric),
         metric = gsub("ETA", "ETAs", metric),
         metric = gsub("RSE", " (RSEs)", metric),
         metric = gsub("B", "b", metric),
         metric = gsub("ka", "Ka", metric),
         metric = factor(metric, levels = unique(metric)),
         parameter = gsub("V", "Vd", parameter),
         parameter = factor(parameter, levels = unique(parameter))) %>% 
  filter(!str_detect(metric, "RSEs"))
pl <- ggplot(dat_long, aes(x = parameter, y = value, colour = parameter)) +
  geom_boxplot() +
  facet_wrap(~metric, scales = "free", nrow = 1) +
  theme_bw() +
  labs(title = "SNP Effects Added (Vd fixed)", x = "Parameter", y = "Value") +
  theme(legend.position = "none",
        strip.text = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, vjust = 0, size = 12, face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
graph2ppt(pl, "Fig_1", append = TRUE, width = 18 * x, height = 9 * x)

# ML methods
heatmap_data <- expand.grid(
  SNPs = c("1,000 SNPs", "10,000 SNPs", "100,000 SNPs", "1,000,000 SNPs"),
  Pruning = c( "LD-based pruning\n+ PCA", "LD-based pruning","No dimensional\nreduction"),
  Method = c("GWAS", "Penalized\nregression", "Random\nForest")
  ) %>%
  mutate(Executed = "Executed",
         Executed = ifelse(SNPs %in% c("100,000 SNPs", "1,000,000 SNPs") & Pruning == "No dimensional\nreduction" & Method != "GWAS", "Not Executed", Executed),
         Executed = ifelse(SNPs == "1,000,000 SNPs" & Method != "GWAS" & Pruning ==  "LD-based pruning", "Not Executed", Executed),
         Executed = ifelse(Method == "GWAS" & Pruning == "LD-based pruning\n+ PCA", "Not Executed", Executed),
         Executed = factor(Executed, levels = c("Executed", "Not Executed")))
pl <- ggplot(heatmap_data, aes(x = Method, y = Pruning, fill = Executed)) +
  geom_tile(color = "white") +
  facet_wrap(~SNPs, nrow = 1) +
  scale_fill_manual(values = c("Executed" = "blue", "Not Executed" = "white")) +
  labs(title = "Covariate Selection Approaches", x = " ", y = " ") +
  theme_minimal() +
  theme(legend.position = "none",
        strip.text = element_text(size = 9),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, vjust = 0, size = 12, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
graph2ppt(pl, "Fig_1", append = TRUE, width = 16 * 0.4, height = 9 * 0.4)

# Imputation performance
dat <- tibble(Method = c("Method A", "Method B", "Method C"),
              `Coverage Rate (% of CIs with reference)` = c(0.91, 0.88, 0.93),  
              `Average CI Width (Statistical Efficiency)` = c(0.2, 0.15, 0.12),
              `Bias (e.g., eMLAR, MRPE, rMPE)` = c(-0.02, 0.01, 0.03),
              `Precision (e.g., eMALAR, RMSRE, MAPE)` = c(0.05, 0.04, 0.03)) %>%
  pivot_longer(cols = -Method, names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = factor(Metric, levels = unique(Metric)))
pl <- ggplot(dat, aes(x = Method, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  geom_hline(data = dat %>% filter(Metric == "Coverage Rate (% of CIs with reference)"), aes(yintercept = 0.9), 
             linetype = "dashed", color = "black", size = 1) +  # 90% Nominal Line for CR
  facet_wrap(~ Metric, scales = "free_y") +
  labs(title = "", x = "", y = "Value") +
  theme_bw() +
  theme(legend.position = "none",
        strip.text = element_text(size = 10),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, vjust = 0, size = 12, face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
graph2ppt(pl, "Fig_1", append = TRUE, width = 16 * 0.4, height = 8 * 0.4)


# GWAS plots
library(ggplotify) # To change manhattan plots to ggplots
manhattan2 <- function (x, chr = "CHR", bp = "BP", p = "P", snp = "SNP", col = c("deepskyblue", "gray10"), suggestiveline = -log10(1e-05), 
                        genomewideline = -log10(5e-08), highlight = NULL, logp = TRUE, 
                        annotatePval = 0.9, annotateTop = FALSE, ...) # annotatePval as a placeholder - we want to annotate the top SNPs
{
  CHR = BP = P = index = NULL
  d = data.frame(CHR = x[[chr]], BP = x[[bp]], P = x[[p]], 
                 pos = NA, index = NA, SNP = x[[snp]], stringsAsFactors = FALSE)
  d <- d[order(d$CHR, d$BP), ]
  d$logp <- -log10(d$P)
  d$index = rep.int(seq_along(unique(d$CHR)), times = tapply(d$SNP, d$CHR, length))
  nchr = length(unique(d$CHR))
  d$pos = d$BP
  xlabel = paste("Chromosome", unique(d$CHR), "position")
  xmax = ceiling(max(d$pos) * 1.03)
  xmin = floor(max(d$pos) * -0.03)
  def_args <- list(xaxt = "n", bty = "n", xaxs = "i", yaxs = "i", 
                   las = 1, pch = 20, xlim = c(xmin, xmax), 
                   ylim = c(0, ceiling(max(d$logp))), xlab = xlabel, ylab = expression(-log[10](italic(p))))
  dotargs <- list(...)
  do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% names(dotargs)]))
  axis(1, ...)
  col = rep_len(col, max(d$index))
  with(d, points(pos, logp, pch = 20, col = col[1], ...))
  if (suggestiveline) abline(h = suggestiveline, col = "blue")
  if (genomewideline) abline(h = genomewideline, col = "red")
  if (!is.null(highlight)) {
    if (any(!(highlight %in% d$SNP))) 
      warning("You're trying to highlight SNPs that don't exist in your results.")
    d.highlight = d[which(d$SNP %in% highlight), ]
    with(d.highlight, points(pos, logp, col = "red", pch = 20, 
                             ...))
  }
  
  
  if (!is.null(annotatePval)) {
    if (logp) {
      topHits = subset(d, SNP %in% paste0("SNP", 1:9))
    }
    else topHits = subset(d, SNP %in% paste0("SNP", 1:9))
    par(xpd = TRUE)
    if (annotateTop == FALSE) {
      if (logp) {
        with(subset(d, SNP %in% paste0("SNP", 1:9)), textxy(pos, 
                                                            -log10(P), offset = 0.625, labs = gsub("SNP", "SNP ", topHits$SNP), 
                                                            cex = 0.7), ...)
      }
      else with(subset(d, SNP %in% paste0("SNP", 1:9)), textxy(pos, 
                                                               P, offset = 0.625, labs = gsub("SNP", "SNP ", topHits$SNP), 
                                                               cex = 0.7), 
                ...)
    }
    else {
      topHits <- topHits[order(topHits$P), ]
      topSNPs <- NULL
      for (i in unique(topHits$CHR)) {
        chrSNPs <- topHits[topHits$CHR == i, ]
        topSNPs <- rbind(topSNPs, chrSNPs[1, ])
      }
      if (logp) {
        textxy(topSNPs$pos, -log10(topSNPs$P), offset = 0.625, 
               labs = gsub("SNP", "SNP ", topSNPs$SNP), cex = 0.5, ...)
      }
      else textxy(topSNPs$pos, topSNPs$P, offset = 0.625, 
                  labs = gsub("SNP", "SNP ", topSNPs$SNP), cex = 0.5, ...)
    }
  }
  par(xpd = FALSE)
}
environment(manhattan2) <- environment(manhattan)
results <- read_csv("high_plot/gwas_10_3_1.csv", show_col_types = FALSE)
pl <- as.ggplot(~ manhattan2(results, highlight = paste0("SNP", 1:9), main = ""))
graph2ppt(pl, "Fig_1", append = TRUE, width = 16 * 0.4, height = 8 * 0.4)


# Cossac-time
time_tb <- read_csv("cossac_test.csv", show_col_types = FALSE) %>%
  mutate(end = ymd_hms(gsub("T ", "", end)),
         start = ymd_hms(gsub("T ", "", start)),
         time_h = as.numeric(end-start, "hours")) %>%
  select(`Number of SNPs` = n_snps, `Time (hours)` = time_h)
pl <- ggplot(time_tb, aes(x = `Number of SNPs`, y = `Time (hours)`)) +
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ poly(x, 3), se = TRUE, color = "blue", linewidth = 0.5) +
  ggtitle("Modeling Runtime of COSSAC-SAMBA by\nNumber of SNPs") +
  theme_bw() + 
  theme(axis.text = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5, vjust = 0, size = 12, face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  annotate("text", x = 400, y = 11,  label = "italic(y) == italic(x)^3", parse = TRUE, colour = "blue", size = 4)
graph2ppt(pl, "Fig_1", append = TRUE, width = 10 * 0.4, height = 10 * 0.4)



