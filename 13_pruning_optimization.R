# Load libraries
library(data.table)
library(tidyverse)
library(ggh4x)
library(ggrepel)

# Load results
pruning_results <- read.table("pruning_results.txt", header = TRUE) %>%
  mutate(N_SNP = gsub("0_", "e", N_SNP),
         P_Retained_SNPs = Retained_SNPs/as.numeric(N_SNP) * 100,
         N_SNP = factor(paste("Total:", format(as.numeric(N_SNP), big.mark = ",", scientific = FALSE), "SNPs")),
         R2 = factor(R2),
         Label = as.character(round(P_Retained_SNPs)),
         Step = paste0("Step size: ", Step, " SNPs"),
         Step = factor(Step, levels = unique(Step)),
         Window = factor(Window, levels = unique(Window))) %>%
  tibble() 

# Compute medians for each Window-R2 group
median_results <- pruning_results %>%
  group_by(N_SNP, Step, Window, R2) %>%
  summarise(median_retained = median(P_Retained_SNPs), .groups = "drop") %>%
  mutate(Label = as.character(round(median_retained)))

# Box plots to visualize SNP retention
pl <- ggplot(pruning_results, aes(x = Window, y = P_Retained_SNPs, color = R2,
                                  group = interaction(Window, R2))) +
  geom_boxplot(aes(fill = R2), alpha = 0.5) +  # Box plots with transparency
  # geom_jitter(width = 0.2, alpha = 0.3) +  # Adds points for better visibility of spread
  geom_text_repel(data = median_results,
                  aes(x = Window, y = median_retained, group = R2, label = Label), 
                  size = 4, box.padding = 0.5, 
                  point.padding = 0.3, max.overlaps = 20, show.legend = FALSE) +  
  geom_line(data = median_results, aes(x = Window, y = median_retained, group = R2, color = R2), 
            size = 1, linetype = "dotted", show.legend = FALSE) +  # Median trend line
  facet_nested_wrap(vars(N_SNP, Step), ncol = 6) +
  theme_bw() +
  labs(title = "Effect of LD Pruning Parameters on SNP Retention (%)",
       x = "Window Size (KB)",
       y = "Retained SNPs (% of Total)",
       color = expression(R^2),
       fill = expression(R^2)) +  # Legend title
  theme(legend.position = "bottom",
        strip.text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0, vjust = 0, size = 14, face = "bold"))
x <- 1
ggsave("prune_optimization.png", height = 9 * x, width = 16 * x)

# Window Size: 1000
# Step-size limited effect: 5 for granularity (slightly slower than 50 but difference is so little)
# R2: 0.1

# The pruned-in and out are already available - not to redo it unless it is fast enough through R.


# No step_size
# Load results
pruning_results <- read.table("pruning_results.txt", header = TRUE) %>%
  filter(Step == 5) %>%
  mutate(N_SNP = gsub("0_", "e", N_SNP),
         P_Retained_SNPs = Retained_SNPs/as.numeric(N_SNP) * 100,
         N_SNP = factor(paste("Total:", format(as.numeric(N_SNP), big.mark = ",", scientific = FALSE), "SNPs")),
         R2 = factor(R2),
         Label = as.character(round(P_Retained_SNPs)),
         Window = factor(Window, levels = unique(Window))) %>%
  tibble() 

# Compute medians for each Window-R2 group
median_results <- pruning_results %>%
  group_by(N_SNP, Window, R2) %>%
  summarise(median_retained = median(P_Retained_SNPs), .groups = "drop") %>%
  mutate(Label = as.character(round(median_retained)))

# Box plots to visualize SNP retention
pl <- ggplot(pruning_results, aes(x = Window, y = P_Retained_SNPs, color = R2,
                                  group = interaction(Window, R2))) +
  geom_boxplot(aes(fill = R2), alpha = 0.5) +  # Box plots with transparency
  # geom_jitter(width = 0.2, alpha = 0.3) +  # Adds points for better visibility of spread
  geom_text_repel(data = median_results,
                  aes(x = Window, y = median_retained, group = R2, label = Label), 
                  size = 4, box.padding = 0.5, 
                  point.padding = 0.3, max.overlaps = 20, show.legend = FALSE) +  
  geom_line(data = median_results, aes(x = Window, y = median_retained, group = R2, color = R2), 
            size = 1, linetype = "dotted", show.legend = FALSE) +  # Median trend line
  facet_wrap(~N_SNP, ncol = 2) +
  theme_bw() +
  labs(title = "Effect of LD Pruning Parameters on SNP Retention (%), with Step Size = 5 SNPs",
       x = "Window Size (KB)",
       y = "Retained SNPs (% of Total)",
       color = expression(R^2),
       fill = expression(R^2)) +  # Legend title
  theme(legend.position = "bottom",
        strip.text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0, vjust = 0, size = 14, face = "bold"))
x <- 0.8
ggsave("prune_optimization_step5.png", height = 9 * x, width = 16 * x)

# Window Size: 1000
# Step-size limited effect: 5 for granularity (slightly slower than 50 but difference is so little)
# R2: 0.1
chosen_params <- "_1000_5_0.1.prune.in"

# Steps: Generate reduced datasets, when including correlated covariates, perform LD analysis using the complete, not reduced datasets
# Generate reduced datasets
n_datasets <- 100
n_snps <- c("10_3", "10_4", "10_5", "10_6")

for (n_snp in n_snps) {
  # Add folders if they don't exist
  folder_path <- paste0("pruned_", n_snp)
  if (!dir.exists(folder_path)) dir.create(folder_path)
  
  for (j in 1:n_datasets) {
    # Only include the subset of pruned SNPs 
    system(paste0("/pub59/iasiimwe/plink1.9/plink --bfile ", n_snp, "/dat_",n_snp, "_", j, 
                  " --extract pruning/pruned_", n_snp, "_", j, chosen_params, " --make-bed --out ",
                  folder_path, "/dat_", n_snp, "_", j))
    message(paste0("NSNPs: ", n_snp, 
                   "; ", round(j*100/n_datasets, 3), "% complete!"))
  }
}
