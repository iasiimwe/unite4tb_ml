# Load required packages
library(tidyverse)
library(ggrepel) 
library(tidytext)
library(ggbump)
library(wesanderson)
library(ggprism)
library(export)
library(ComplexUpset)
library(patchwork)
library(RColorBrewer)

# Correlated function
gsub_x <- function(x) gsub("SNP\\d, ", "", x)
methods_to_plot <- c("glmnet", "rf", "gwas", "gwas_glmnet", "gwas_rf")

# Formatting funcion
format_for_plot <- function(dat) {
  dat %>%
    arrange(desc(F1_score)) %>%
    mutate(Method = case_when(Method == "gwas" ~ "GWAS",
                              Method == "rf" ~ "Random Forest",
                              Method == "glmnet" ~ "Penalized\nregression",
                              Method == "gwas_glmnet" ~ "GWAS + Penalized\nregression",
                              Method == "gwas_rf" ~ "GWAS + Random\nForest"),
           F1_score = round(F1_score, 3),
           Method2 = paste0(Method, "\n(", str_pad(as.character(F1_score), width = 5, pad = "0", "right"), ")"),
           Method = fct_rev(factor(Method, levels = unique(Method)))) %>%
    select(-F1_score) %>%
    separate(Covars, into = paste0("SNP", 1:9), sep = "; ") %>%
    mutate_at(vars(SNP1:SNP9), gsub_x) %>%
    mutate_at(vars(SNP1:SNP9), as.numeric) %>%
    pivot_longer(SNP1:SNP9, names_to = "Covariates", values_to = "Percent selection") %>%
    arrange(Method, desc(`Percent selection`)) %>%
    mutate(Covariates = factor(Covariates, levels = unique(Covariates)),
           label = as.character(`Percent selection`),
           methodnumeric = as.numeric(Method)) %>%
    return()
}

# Bumpplot function
bumpplot_fn <- function(dat, ptitle = "A. Low Effect Size") {
  # display.brewer.all()
  green_colours <- brewer.pal(5, "Greens")
  blue_colours <- brewer.pal(5, "Blues")
  red_colours <- brewer.pal(5, "Reds")
  manual_colours <- c(SNP9 = green_colours[5], SNP8 = green_colours[3], SNP7 = green_colours[2],
                      SNP6 = blue_colours[5], SNP5 = blue_colours[3], SNP4 = blue_colours[2],
                      SNP3 = red_colours[5], SNP2 = red_colours[3], SNP1 = red_colours[2])
  ggplot(dat, aes(x = methodnumeric , y = `Percent selection`)) +
    geom_bump(aes(color = Covariates, group=Covariates), linewidth = 1, alpha = 0.9, show.legend = FALSE) +
    geom_point(aes(color = Covariates), size = 1, alpha = 0.9, show.legend = FALSE) +
    geom_label_repel(aes(label = label, color = Covariates), size = 3, 
                     fontface = "bold", alpha = 0.9, direction = "y", max.overlaps = 40) + 
    geom_label_repel(data = dat %>% filter(methodnumeric == 1), # filter(methodnumeric == length(unique(results$Method))
                     aes(label = Covariates, color = Covariates),
                     size = 3, hjust = 1, show.legend = FALSE, direction = "y", nudge_x = -0.2) +
    facet_wrap(~Cor, ncol = 1) +
    scale_size_area(max_size = 20) +
    scale_alpha_continuous(range = c(0.4, 0.9)) + 
    guides(size = "none", color = "none", alpha = "none", shape = "none") +
    theme_bw() +
    theme(legend.position = "null",
          strip.text = element_text(size = 12, face = "bold"),
          axis.title.y = element_text(size = 12, face = "bold"),
          axis.title.x = element_text(size = 12, face = "bold"),
          axis.ticks.x = element_blank(), 
          axis.text = element_text(size = 10, colour = "black"),
          plot.title = element_text(hjust = 0, vjust = 0, size = 13, face = "bold")) +
    labs(x = "Methods (F1 scores based on correlations)", 
         y = "Percent covariate selection",
         title = ptitle) +
    scale_x_continuous(
      limits = c(0.8, max(dat$methodnumeric) + 0.05),
      breaks = unique(dat$methodnumeric),
      labels = unique(dat$Method2),
      guide = guide_axis(n.dodge = 1)) +
    coord_cartesian(ylim = c(30, 100)) +
    scale_color_manual(values = manual_colours) %>%
    return()
}

# Other relevant functions
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

summary_not_correlated_1000 <- function(path) {
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

summary_correlated_1000 <- function(path) {
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

# Performance
# -----------
pl_fn <- function(df) {
  df %>%
    mutate(SD = parse_number(gsub(".*\\(([^)]+)\\).*", "\\1", F1)),
           F1 = as.numeric(gsub(" .*", "", F1)),
           upper = F1 + qnorm(0.975) * SD/sqrt(n_datasets),
           lower = F1 + qnorm(0.025) * SD/sqrt(n_datasets),
           Scenario = paste0(str_to_title(Scenario), " effect size"),
           Scenario = factor(Scenario, levels = c("Low effect size", "High effect size")),
           NSNPs = gsub("0_", "e", NSNPs),
           NSNPs = paste(format(as.numeric(NSNPs), big.mark = ",", scientific = FALSE), "SNPs"),
           Pruned = case_when(Pruned == "none" ~ "None", Pruned == "pruned" ~ "Pruned", 
                              Pruned == "pca" ~ "Pruned + PCA", Pruned == "pca_o" ~ "PCA"),
           Pruned = factor(Pruned, levels = c("None", "Pruned", "Pruned + PCA", "PCA")),
           label = as.character(round(F1, 2)),
           label = if_else(str_detect(label, "\\."), label, paste0(label, ".00")),
           label = if_else(nchar(label) == 4, label, paste0(label, "0")),
           label = if_else(perfect == 0, label, paste0(label, " (", perfect, ")")),
           perfect = if_else(perfect == 0, NA, perfect),
           NSNPs = factor(NSNPs),
           Method = case_when(Method == "gwas" ~ "GWAS",
                              Method == "rf" ~ "Random Forest",
                              Method == "glmnet" ~ "Penalized\nregression",
                              Method == "gwas_glmnet" ~ "GWAS + Penalized\nregression",
                              Method == "gwas_rf" ~ "GWAS + Random\nForest"),
           Method = factor(Method, levels = c("GWAS", "Random Forest", "GWAS + Random\nForest", "GWAS + Penalized\nregression", "Penalized\nregression"))) %>%
    ggplot(aes(x = Method, y = F1, color = Scenario)) +
    geom_point(size = 1.5, position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, position = position_dodge(width = 0.5)) +
    geom_text(aes(label = label, hjust = ifelse(Scenario == "Low effect size", 1.2, -0.2),
                  vjust = ifelse(Scenario == "Low effect size", 1.7, -0.5)), 
              size = 3, show.legend = FALSE) +
    facet_grid(Pruned ~ NSNPs) +
    labs(title = "A. All SNP sizes (correlated covariates not accounted for)",
         x = "Method",
         y = "F1 Score",
         color = "Scenario") +
    theme_bw() + 
    theme(legend.position = "bottom",
          strip.text = element_text(size = 12, face = "bold"),
          axis.text = element_text(size = 10, colour = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title = element_text(size = 12, face = "bold"),
          legend.title = element_text(size = 12, face = "bold"),
          plot.title = element_text(hjust = 0, vjust = 0, size = 12, face = "bold"),
          # panel.grid.major = element_blank(), 
          # panel.grid.minor = element_blank()
          ) +
    coord_cartesian(ylim = c(0, 1)) %>%
    return()
}
pl_fn_1000 <- function(df) {
  df %>%
    mutate(SD = parse_number(gsub(".*\\(([^)]+)\\).*", "\\1", F1)),
           F1 = as.numeric(gsub(" .*", "", F1)),
           upper = F1 + qnorm(0.975) * SD/sqrt(n_datasets),
           lower = F1 + qnorm(0.025) * SD/sqrt(n_datasets),
           Scenario = paste0(str_to_title(Scenario), " effect size"),
           Scenario = factor(Scenario, levels = c("Low effect size", "High effect size")),
           NSNPs = gsub("0_", "e", NSNPs),
           NSNPs = paste(format(as.numeric(NSNPs), big.mark = ",", scientific = FALSE), "SNPs"),
           label = as.character(round(F1, 2)),
           label = if_else(str_detect(label, "\\."), label, paste0(label, ".00")),
           label = if_else(nchar(label) == 4, label, paste0(label, "0")),
           label = if_else(perfect == 0, label, paste0(label, " (", perfect, ")")),
           perfect = if_else(perfect == 0, NA, perfect),
           NSNPs = factor(NSNPs),
           Cor = factor(Cor, levels = unique(dat$Cor)),
           Method = case_when(Method == "gwas" ~ "GWAS",
                              Method == "rf" ~ "Random Forest",
                              Method == "glmnet" ~ "Penalized\nregression",
                              Method == "gwas_glmnet" ~ "GWAS + Penalized\nregression",
                              Method == "gwas_rf" ~ "GWAS + Random\nForest"),
           Method = factor(Method, levels = c("GWAS", "Random Forest", "GWAS + Random\nForest", "GWAS + Penalized\nregression", "Penalized\nregression"))) %>%
    ggplot(aes(x = Method, y = F1, color = Scenario)) +
    geom_point(size = 1.5, position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, position = position_dodge(width = 0.5)) +
    geom_text(aes(label = label, hjust = ifelse(Scenario == "Low effect size", 1.2, -0.2),
                  vjust = ifelse(Scenario == "Low effect size", 1.7, -0.5)), 
              size = 3, show.legend = FALSE) +
    facet_wrap(~Cor, nrow = 1) +
    labs(title = "B. N = 1000 SNPs",
         x = "Method",
         y = "F1 Score",
         color = "Scenario") +
    theme_bw() + 
    theme(legend.position = "bottom",
          strip.text = element_text(size = 12, face = "bold"),
          axis.text = element_text(size = 10, colour = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title = element_text(size = 12, face = "bold"),
          legend.title = element_text(size = 12, face = "bold"),
          plot.title = element_text(hjust = 0, vjust = 0, size = 12, face = "bold"),
          # panel.grid.major = element_blank(), 
          # panel.grid.minor = element_blank()
          ) +
    coord_cartesian(ylim = c(0, 1)) %>%
    return()
}

n_datasets <- 100
n_round <- 3 

# Not_correlated (all)
dat <- summary_not_correlated("Eta_sim_performance1.csv") 
pl <- pl_fn(dat) # + coord_cartesian(ylim = c(0, 1.1))
x <- 0.8
graph2ppt(pl, "Fig4", append = FALSE, width = 16 * x, height = 9.5 * x)

# Not_correlated + correlated (1000)
dat <- mutate(summary_not_correlated_1000("Eta_sim_performance1_1000.csv"), Cor = "Correlated covariates not accounted for") %>%
  bind_rows(mutate(summary_correlated_1000("Eta_sim_performance2_1000.csv"), Cor = "Correlated covariates accounted for (R² > 0.5)")) %>%
  bind_rows(mutate(summary_correlated_1000("Eta_sim_performance1_1000.csv"), Cor = "Correlated covariates accounted for (R² > 0.1)"))
pl <- pl_fn_1000(dat) 
graph2ppt(pl, "Fig4", append = TRUE, width = 16 * x, height = 5 * x)


## Figure 5
# ---------
datA <- read_csv("Eta_sim_performance1_1000.csv", show_col_types = FALSE) %>% 
  select(Method, F1_score = `F1 score mean`, Covars = `True covariate % selection`, Perfect = `Perfect selection %`, Scenario) %>%
  format_for_plot() %>%
  mutate(Cor = "Correlated covariates not accounted for")
datB <- read_csv("Eta_sim_performance2_1000.csv", show_col_types = FALSE) %>% 
  select(Method, F1_score = `F1 score cor mean`, Covars = `Correlated covariate % selection`, Perfect = `Perfect selectioncor %`, Scenario) %>%
  format_for_plot() %>%
  mutate(Cor = "Correlated covariates accounted for (R² = 0.5)")
datC <- read_csv("Eta_sim_performance1_1000.csv", show_col_types = FALSE) %>% 
  select(Method, F1_score = `F1 score cor mean`, Covars = `Correlated covariate % selection`, Perfect = `Perfect selectioncor %`, Scenario) %>%
  format_for_plot() %>%
  mutate(Cor = "Correlated covariates accounted for (R² = 0.1)")
dat <- datA %>%
  bind_rows(datB) %>%
  bind_rows(datC)
method2_dat <- dat %>%
  mutate(Method2 = gsub("[()]", "", str_extract(Method2, "\\(.*\\)"))) %>%
  group_by(Method, Scenario) %>%
  summarise(Method2 = paste0(sort(unique(paste0(Method2, gsub("[0-9.]", "", Cor)))), collapse = ","), .groups = "drop") %>%
  mutate(Method2 = gsub("[^0-9.,]", "", Method2),
         Method2 = gsub(",", ", ", Method2),
         Method2 = paste0(Method, "\n(", Method2, ")"))
dat <- dat %>%
  select(-Method2) %>%
  left_join(method2_dat) %>%
  mutate(Cor = factor(Cor, levels = c("Correlated covariates not accounted for",
                                      "Correlated covariates accounted for (R² > 0.5)",
                                      "Correlated covariates accounted for (R² > 0.1)")))

scenarios <- c("low", "high")
for (i in seq_along(scenarios)) {
  plot_title <- ifelse(scenarios[[i]] == "low", "A. Low Effect Size", "B. High Effect Size")
  bumpplot <- dat %>%
    filter(Scenario == scenarios[[i]]) %>%
    bumpplot_fn(plot_title)
  x <- 0.8
  append_status <- ifelse(i == 1, FALSE, TRUE)
  graph2ppt(bumpplot, file = "Fig5", append = append_status, width = 16 * x, height = 12 * x) 
}


# Legend
# Define SNPs and corresponding attributes
SNPs <- c("SNP1", "SNP2", "SNP3", "SNP4", "SNP5", "SNP6", "SNP7", "SNP8", "SNP9")
MAF <- rep(c("5%", "10%", "20%"), times = 3)  # Repeated for each column
Correlation <- rep(c("Low", "Moderate", "High"), each = 3)  # Repeated for each row

# Define colors
green_colours <- brewer.pal(5, "Greens")
blue_colours <- brewer.pal(5, "Blues")
red_colours <- brewer.pal(5, "Reds")

manual_colours <- c(SNP9 = green_colours[5], SNP8 = green_colours[3], SNP7 = green_colours[2],
                    SNP6 = blue_colours[5], SNP5 = blue_colours[3], SNP4 = blue_colours[2],
                    SNP3 = red_colours[5], SNP2 = red_colours[3], SNP1 = red_colours[2])

# Create a dataframe
df <- tibble(SNPs, MAF, Correlation) %>%
  mutate(MAF = factor(MAF, levels = unique(MAF)),
         Correlation = factor(Correlation, levels = unique(Correlation)))

# Heatmap plot
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
x <- 0.3
graph2ppt(pl, file = "Fig5_legend", append = FALSE, width = 12 * x, height = 10 * x) 


## Time
# -----
str_pad_fn <- function(x, n_round = 1) {
  x <- as.character(round(x, n_round))
  x <- ifelse(!str_detect(x, "\\."), paste0(x, "."), x)
  
  x <- str_pad(x, width = n_round + nchar(gsub("\\..*", "", x)) + 1, pad = "0", "right")
  return(x)
} 
dat <- read_csv("Time.csv", show_col_types = FALSE) %>%
  mutate(mean = mean/60,
         sd = sd/60,
         upper = mean + qnorm(0.975) * sd/sqrt(n_datasets),
         lower = mean + qnorm(0.025) * sd/sqrt(n_datasets),
         Scenario = paste0(str_to_title(Scenario), " effect size"),
         NSNPs = gsub("0_", "e", NSNPs),
         NSNPs = paste(format(as.numeric(NSNPs), big.mark = ",", scientific = FALSE), "SNPs"),
         NSNPs = factor(NSNPs),
         Pruned = case_when(Pruned == "none" ~ "None", Pruned == "pruned" ~ "Pruned", 
                            Pruned == "pca" ~ "Pruned + PCA", Pruned == "pca_o" ~ "PCA"),
         Pruned = factor(Pruned, levels = c("None", "Pruned", "Pruned + PCA", "PCA")),
         Method = case_when(Method == "gwas" ~ "GWAS",
                            Method == "rf" ~ "Random Forest",
                            Method == "glmnet" ~ "Penalized\nregression",
                            Method == "gwas_glmnet" ~ "GWAS + Penalized\nregression",
                            Method == "gwas_rf" ~ "GWAS + Random\nForest"),
         Method = factor(Method, levels = c("GWAS", "Random Forest", "GWAS + Random\nForest", "GWAS + Penalized\nregression", "Penalized\nregression")),
         label = if_else(mean < 1, 
                         paste0(str_pad_fn(mean * 60), "±\n", str_pad_fn(sd * 60), "s"),
                         paste0(str_pad_fn(mean), "±\n", str_pad_fn(sd), "m"))) 
pl <- dat %>%
  ggplot(aes(x = Method, y = mean, color = Scenario)) +
  geom_point(size = 2.5, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, position = position_dodge(width = 0.5)) +
  geom_text(aes(label = label, hjust = ifelse(Scenario == "High effect size", 1.2, -0.2),
                vjust = ifelse(mean < 10, -0.3, 1.5)), 
            size = 3.5, show.legend = FALSE) +
  facet_grid(Pruned ~ NSNPs) +
  labs(title = "Computational cost/time",
       x = "Method",
       y = "Time",
       color = "Scenario") +
  theme_bw() + 
  theme(legend.position = "bottom",
        strip.text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10, colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0, vjust = 0, size = 12, face = "bold")) 
x <- 0.9
graph2ppt(pl, file = "ML_cov_time", append = FALSE, width = 19 * x, height = 9 * x) 







