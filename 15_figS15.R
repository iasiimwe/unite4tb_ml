# Load required packages
# ---------------------
library(tidyverse)

# Single vs Multiple Imputation
# Rubins rules
rubin_fn <- function(coef_estimates, coef_ses) {
  m <- length(coef_estimates)
  qbar <- mean(coef_estimates)
  ubar <- mean(coef_ses^2)
  B <- var(coef_estimates)
  pooled_var <- ubar + (1 + (1/m)) * B
  pooled_se <- sqrt(pooled_var)
  df <- (m - 1) * ((ubar + (1 + 1/m) * B)/((1 + 1/m) * B))^2
  t_stat <- qbar / pooled_se
  
  # Safe log p-value calculation
  log_p <- pt(abs(t_stat), df, lower.tail = FALSE, log.p = TRUE)
  p_val <- 2 * exp(log_p)
  
  # Format for display only
  p_val_fmt <- ifelse(p_val < .Machine$double.xmin,
                      format(.Machine$double.xmin, scientific = TRUE),
                      format(p_val, scientific = TRUE))
  
  return(list(
    pooled_est = qbar,
    pooled_se = pooled_se,
    df = df,
    t = t_stat,
    p = p_val,        # numeric
    p_fmt = p_val_fmt # character (for display)
  ))
}

# Artifical dataset
set.seed(7)  # for reproducibility

# Simulate genotype as 0 or 1 (e.g., dominant coding)
Best_Guess_Genotype <- rbinom(100, 1, prob = 0.5)

# Simulate outcome with some noise
Outcome <- 3.52 + 1.16 * Best_Guess_Genotype + rnorm(100, mean = 0, sd = 1) # apparent oral clearance (CL) = 3.52, low effect size

# Create subject IDs
ID <- 1:100
dat <- tibble(ID, Outcome, Best_Guess_Genotype)

# Fit the model
model <- lm(Outcome ~ Best_Guess_Genotype, data = dat)
p_value <- summary(model)$coefficients["Best_Guess_Genotype", "Pr(>|t|)"] # 2.306349e-09
beta <- summary(model)$coefficients["Best_Guess_Genotype", "Estimate"] # 1.173997
se <- summary(model)$coefficients["Best_Guess_Genotype", "Std. Error"] # 0.1783729


# Assuming a random 10% of patients with high outcomes and genotype 1, had probabilities of (0.1, 0.9, 0)
# I.e. in 10% of the datasets, their probability would be 0, not 1.
set.seed(7)
ids <- dat %>%
  filter(Outcome > median(dat$Outcome), Best_Guess_Genotype == 1) %>%
  slice_sample(n = 10) %>%
  pull(ID)

imp <- 10 # Number of 'imputed' datasets
betas <- rep(beta, imp)
ses <- rep(se, imp)

dat_10_percent <- dat 
dat_10_percent$Best_Guess_Genotype[ids] <- 0
model2 <- lm(Outcome ~ Best_Guess_Genotype, data = dat_10_percent)
betas[1] <- summary(model2)$coefficients["Best_Guess_Genotype", "Estimate"] # 0.8265813
ses[1] <- summary(model2)$coefficients["Best_Guess_Genotype", "Std. Error"] # 0.2050813
rubin_fn(betas, ses)$p # 6.000419e-07   # Single imputation would over-estimate
rubin_fn(betas, ses)$pooled_est # 1.139256
rubin_fn(betas, ses)$pooled_se # 0.2147506


# Assuming a random 10% of patients with high outcomes and genotype 0, had probabilities of (0.9, 0.1, 0)
# I.e. in 10% of the datasets, their probability would be 1, not 0.
set.seed(7)
ids <- dat %>%
  filter(Outcome > median(dat$Outcome), Best_Guess_Genotype == 0) %>%
  slice_sample(n = 10) %>%
  pull(ID)

imp <- 10 # Number of 'imputed' datasets
betas <- rep(beta, imp)
ses <- rep(se, imp)

dat_10_percent <- dat 
dat_10_percent$Best_Guess_Genotype[ids] <- 1
model2 <- lm(Outcome ~ Best_Guess_Genotype, data = dat_10_percent)
betas[1] <- summary(model2)$coefficients["Best_Guess_Genotype", "Estimate"] # 1.324087
ses[1] <- summary(model2)$coefficients["Best_Guess_Genotype", "Std. Error"] # 0.1695281
rubin_fn(betas, ses)$p # 1.460218e-10
rubin_fn(betas, ses)$pooled_est # 1.189006
rubin_fn(betas, ses)$pooled_se # 0.184356


# Create summary data
plot_dat <- tibble(Method = c("Single", "Multiple (Single Underestimating)", 
                              "Multiple (Single Overestimating)"),
                   p_value = -log10(c(2.306349e-09, 1.460218e-10, 6.000419e-07)),
                   Beta = c(1.173997, 1.189006, 1.139256),
                   SE = c(0.1783729, 0.184356, 0.2147506)) %>%
  mutate_at(vars(Beta, SE), ~round(., 2)) %>%
  mutate(label = paste0("β = ", Beta, " ± ", SE),
         Method = factor(Method, levels = unique(Method)))

ggplot(plot_dat, aes(x = Method, y = p_value, fill = Method)) +
  geom_col(position = "dodge") +
  geom_text(aes(label = label), 
            vjust = -0.5, size = 4, fontface = "italic") +
  labs(title = "Single vs Multiple Imputation P-values",
       y = "-log10(P)",
       x = "Imputation Method",
       fill = "Imputation Method") +
  geom_hline(yintercept = -log10(5e-8), linetype = "dashed", color = "red") +
  theme_bw() +
  theme(axis.text = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0, vjust = 0, size = 12, face = "bold")) +
  coord_cartesian(ylim = c(0, 10.5))
x <- 0.5
ggsave("Single_vs_multiple.png", height = 9 * x, width = 16 * x)


