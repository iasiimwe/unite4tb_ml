# Use our paper: https://ascpt.onlinelibrary.wiley.com/doi/10.1002/cpt.3229
# ------------
# cd /pub59/iasiimwe/TB/ukb; source /pub59/iasiimwe/miniconda3/bin/activate base

# Serum urate
# Chromosome 4, which had SLC2A9, top SNP rs938564
# Try QTL and vQTL analyses and see if results will be different - hopefully by reducing false-positives you reduce the testing burden

# Start with a random 5000 and try to reproduce any key SNPs

# Standard (chromosome 4, 5000 random participants)
# cd /pub59/iasiimwe/TB/ukb
## Obtain genotype data in PLINK format
# cohort="dis"
# outcome="urate"
# # Randomly sample 11000 participants (plan is for 5000, some have no outcome)
# shuf --random-source=<(yes 42) -n 11000 /pub59/iasiimwe/UKBB/GWAS/keep_BFZ_${cohort}.txt > ukb_sampled_participants.txt


# for i in 4 6; do
# # Convert to vcf
#  /pub59/iasiimwe/plink2 \
#   --bgen /pub59/iasiimwe/ukbb_geno/ukb22828_c${i}_b0_v3.bgen ref-first \
#   --sample /pub59/iasiimwe/ukbb_geno/ukb22828_c${i}_b0_v3_s487159.sample \
#   --keep ukb_sampled_participants.txt \
#   --export vcf bgz vcf-dosage=GP \
#   --out ${cohort}_chr_${i} \
#   --geno 0.05 \
#   --maf 0.05 \
#   --hwe 0.00001 \
#   --rm-dup 'force-first'
#
# # Hard genotypes
#  /pub59/iasiimwe/plink2 \
#   --bgen /pub59/iasiimwe/ukbb_geno/ukb22828_c${i}_b0_v3.bgen ref-first \
#   --sample /pub59/iasiimwe/ukbb_geno/ukb22828_c${i}_b0_v3_s487159.sample \
#   --keep ukb_sampled_participants.txt \
#   --make-bed \
#   --out ${cohort}_chr_${i} \
#   --geno 0.05 \
#   --maf 0.05 \
#   --hwe 0.00001 \
#   --minimac3-r2-filter 0.3 \
#   --hard-call-threshold 0.1 \
#   --rm-dup 'force-first'
# done
#
# # Convert VCF v4.3 to v4.2 and get probabilities
# for i in 4 6; do
# zcat ${cohort}_chr_${i}.vcf.gz | sed 's/VCFv4.3/VCFv4.2/' | bgzip > converted_${cohort}_chr_${i}.vcf.gz
# /pub59/iasiimwe/tabix/tabix-0.2.6/tabix -p vcf converted_${cohort}_chr_${i}.vcf.gz
#
# # Get probabilities
# /pub59/iasiimwe/vcftools_0.1.13/bin/vcftools \
#   --gzvcf converted_${cohort}_chr_${i}.vcf.gz \
#   --out ${cohort}_chr_${i} \
#   --extract-FORMAT-info GP
#
# # Get genotypes (only imputed have probabilities)
# /pub59/iasiimwe/vcftools_0.1.13/bin/vcftools \
#   --gzvcf converted_${cohort}_chr_${i}.vcf.gz \
#   --out ${cohort}_chr_${i} \
#   --extract-FORMAT-info GT
# done
#
#
# # Get previously-reported SNPs (Table S5 of our paper)
# Base-pair positions for SLC2A9: 827848..10041894 (Assembly: GRCh37.p13 (GCF_000001405.25), https://www.ncbi.nlm.nih.gov/gene/56606)
# BP for SLC17A1: 25783143..25832280 (Assembly: GRCh37.p13 (GCF_000001405.25), https://www.ncbi.nlm.nih.gov/gene/6568)
# Add 10,000 BPs on either side i.e. 
  # SLC2A9: 817848..10051894
  # SLC17A1: 25773143..25842280
#


# # Get plink files
# for i in 4 6; do
# if [ "$i" -eq 4 ]; then bp1=817848; bp2=10051894; else bp1=25773143; bp2=25842280; fi
# /pub59/iasiimwe/plink2 \
#  --bfile ${cohort}_chr_${i} \
#  --chr ${i} --from-bp ${bp1} --to-bp ${bp2} \
#  --geno 0.05 --maf 0.01 --hwe 0.000001 --make-bed \
#  --out ukb_chr_${i}
# done


# Get phenotype and covariates (R) - Focus on urate
# -----------------
library(tidyverse)
library(data.table)
cohort <- "dis"
# Pheno
dat <- fread(paste0("/pub59/iasiimwe/UKBB/GWAS/pheno_BFZ_", cohort, ".txt")) %>%
  filter(FID %in% fread("ukb_sampled_participants.txt")$V1) %>%
  select(FID, IID, urate) # Focus on urate for now
dat <- dat %>%
  filter(!is.na(urate)) 

# Covariates
covar_dat <- fread(paste0("/pub59/iasiimwe/UKBB/GWAS/covar_BFZ_", cohort, ".txt")) %>%
  filter(FID %in% fread("ukb_sampled_participants.txt")$V1) %>%
  select(FID, IID, age, sex, C1:C10) # https://ascpt.onlinelibrary.wiley.com/doi/10.1002/cpt.3229


nsamples <- c(1000, 500, 100)
set.seed(7)
for (nsample in nsamples) {
  dat <- dat %>%
    slice_sample(n = nsample)
  write.table(dat, paste0("ukb_", nsample, "_pheno.txt"), sep = " ", 
              row.names = FALSE, quote = FALSE)
  
  write.table(filter(covar_dat, FID %in% dat$FID), paste0("ukb_", nsample, "_covar.txt"), sep = " ", 
              row.names = FALSE, quote = FALSE)
}


# # Get Info scores
# wget  -nd  biobank.ndph.ox.ac.uk/ukb/ukb/auxdata/ukb_imp_mfi.tgz
# tar -xvzf ukb_imp_mfi.tgz


# Subset GT and GP files in R
# ---------------------------
library(tidyverse)
library(data.table)

n_sample <- 1000
for (i in c(4)) {
  # Read bim file to get the applicable SNPs
  bim <- fread(paste0("ukb_chr_", i, ".bim")) %>%
    select(SNP = V2, POS = V4)
  
  # Sampled participants
  ids <- fread(paste0("ukb_", n_sample, "_pheno.txt"))$FID
  
  # Get info scores
  info <- fread(paste0("ukb_mfi_chr", i, "_v3.txt"))
  info <- info[, .(SNP = V2, POS = V3, R2 = V8)]
  info <- info[SNP %in% bim$SNP]
  info <- info[, .SD[1], by = SNP]
  # nrow(info) # 24509
  fwrite(info, paste0("chr", i, ".info.txt"))
  
  # Genotypes
  imp_dat <- fread(paste0("dis_chr_", i, ".GT.FORMAT")) %>%
    select(-CHROM) %>%
    filter(POS %in% bim$POS)
  colnames(imp_dat) <- gsub("_.*", "", colnames(imp_dat))
  imp_dat <- imp_dat %>% select(any_of(c("POS", ids)))
  fwrite(imp_dat, paste0("ukb_chr_", i, ".GT.txt"))
  
  # Genotype probabilities
  imp_dat <- fread(paste0("dis_chr_", i, ".GP.FORMAT")) %>%
    select(-CHROM) %>%
    filter(POS %in% bim$POS)
  colnames(imp_dat) <- gsub("_.*", "", colnames(imp_dat))
  imp_dat <- imp_dat %>% select(any_of(c("POS", ids)))
  fwrite(imp_dat, paste0("ukb_chr_", i, ".GP.txt"))
}


# Write the datasets then reproduce the above - you will have the betas and SEs
# -------------------------------------------
# Do for 100, but also later test the first 10
# Do for 1000 patients and subset later
library(tidyverse)
library(data.table)

# nsample <- 1000

# Add folders if they don't exist
if (!dir.exists("imputed_datasets")) dir.create("imputed_datasets")

# Convert genotype function
convert_genotype <- function(x) c("A A", "A T", "T T")[x + 1]

# Additive coding function 
additive_coding_fn <- function(geno, r2, pmissing) {
  # Step 1: Split genotype strings and convert to matrix
  prob_mat <- do.call(rbind, strsplit(geno, ",")) |> apply(2, as.numeric)
  
  # Step 2: Compute uncertainty factor
  uncertainty_factor <- (1 - r2) * pmissing
  
  # Step 3: Add noise per row (simulate imputation uncertainty)
  noisy_probs <- prob_mat + matrix(rnorm(length(geno) * 3, mean = 0, sd = rep(uncertainty_factor, each = 3)), 
                                   ncol = 3, byrow = FALSE)
  
  # Step 4: Min-max normalize each row
  row_min <- apply(noisy_probs, 1, min)
  row_max <- apply(noisy_probs, 1, max)
  norm_probs <- (noisy_probs - row_min) / (row_max - row_min + 1e-8)  # add small value to avoid 0/0
  norm_probs <- norm_probs / rowSums(norm_probs)
  
  # Step 5: Sample genotype (0,1,2) from noisy probs
  sampled_genotypes <- apply(norm_probs, 1, function(p) sample(0:2, 1, prob = p))
  
  # Step 6: Convert to PED-style allele pairs (A A, A T, T T)
  return(convert_genotype(sampled_genotypes))
} 

# Do per chromosome
chrs <- c(4)
nsim <- 100
for (chr in chrs) {
  # Get the genotyped data and reshape from wide to long (gather samples under an "ID" column)
  geno_dat <- fread(paste0("ukb_chr_", chr, ".GT.txt"))
  geno_dat <- melt(geno_dat, id.vars = "POS", variable.name = "ID", value.name = "Genotype")
  geno_dat[, Genotype := fifelse(Genotype == "./.", NA_integer_,
                                 fifelse(Genotype == "0/0", 0L,
                                         fifelse(Genotype == "0/1", 1L,
                                                 fifelse(Genotype == "1/1", 2L, NA_integer_))))]
  # Add info score
  info <- fread(paste0("chr", chr, ".info.txt")) 
  geno_dat <- merge(geno_dat, info, by = "POS", all.x = TRUE)

  # Retain only  ID, Genotype and SNP
  geno_dat <- geno_dat %>%
    distinct(ID, SNP, Genotype)
  
  # Get the imputed data and reshape from wide to long
  imp_dat <- fread(paste0("ukb_chr_", chr, ".GP.txt")) %>%
    melt(id.vars = "POS", variable.name = "ID", value.name = "Genotype") %>%
    filter(Genotype != ".")
  imp_dat <- merge(imp_dat, info, by = "POS", all.x = TRUE)

  # Get map files
  imp_map <- imp_dat %>%
    distinct(SNP, POS) %>%
    rename(BP = POS) %>%
    mutate(CHR = chr, .before = SNP) %>%
    mutate(Dist = 0, .before = BP) %>% 
    filter(!is.na(SNP))
  
  # Get ped files
  ped_start <- imp_dat %>%
    distinct(ID) %>% # Placeholder for FID
    mutate(IID = ID, PID = 0, MID = 0, Sex = 0, Pheno = -9)
  
  imp_dat <- imp_dat %>%
    select(ID, Genotype, SNP, R2, POS)
  
  # Get the original genotypes (we can use this earlier to make pmissing (=1) larger for imputed)
  original_snps <- fread(paste0("/pub59/iasiimwe/ukbb_geno/ukb_snp_chr", chr, "_v2.bim"))$V2
  
  set.seed(7)
  for (i in 1:nsim) {
    
    dat_all <- ped_start
    for (j in seq_along(imp_map$SNP)) {
      snp <- imp_map$SNP[[j]]
      geno_dat_snp <- geno_dat[SNP == snp]
      dat_snp <- imp_dat[SNP == snp]
      
      # Missingness percentage
      pmissing <-  ifelse(snp %in% original_snps, 
                          (length(dat_snp$ID) + geno_dat_snp %>% filter(is.na(Genotype)) %>% nrow()) / length(geno_dat_snp$ID),
                          1)

      # Add coded genotypes
      dat_snp[, Coded := additive_coding_fn(Genotype, R2, pmissing)]
      dat_snp <- dat_snp[, .SD[1], by = ID]
      
      # Genotype data
      geno_dat_snp[, Genotype := convert_genotype(Genotype)]
      geno_dat_snp[is.na(Genotype), Genotype := "0 0"]
      geno_dat_snp <- geno_dat_snp[, .SD[1], by = ID]
      
      # Add to dat_all
      dat_all <- merge(dat_all, dat_snp[, .(ID, Coded)], by = "ID", all.x = TRUE) # Imputed data
      set(dat_all, which(is.na(dat_all$Coded)), "Coded", "0 0")
      dat_all <- merge(dat_all, geno_dat_snp[, .(ID, Genotype)], by = "ID", all.x = TRUE) # Genotyped data
      dat_all[, Coded := fifelse(Coded == "0 0", Genotype, Coded)]
      dat_all[, Genotype := NULL]
      setnames(dat_all, "Coded", snp)
      
      message(paste0("Chr: ", chr, "; dataset ", i, " of ", nsim, "; ", round(j * 100/length(imp_map$SNP), 3), "% complete!"))
    }
    
    # Save map and ped files
    dat_all <- dat_all[, .SD[1], by = ID]
    write.table(imp_map, file = paste0("imputed_datasets/chr_", chr, "_dat_", i, ".map"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
    write.table(dat_all, file = paste0("imputed_datasets/chr_", chr, "_dat_", i, ".ped"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  }
}


# Write missingness percentage
library(tidyverse)
library(data.table)

# nsample <- 1000

# Add folders if they don't exist
if (!dir.exists("imputed_datasets")) dir.create("imputed_datasets")

# Do per chromosome
chrs <- c(4)
for (chr in chrs) {
  # Get the genotyped data and reshape from wide to long (gather samples under an "ID" column)
  geno_dat <- fread(paste0("ukb_chr_", chr, ".GT.txt"))
  geno_dat <- melt(geno_dat, id.vars = "POS", variable.name = "ID", value.name = "Genotype")
  geno_dat[, Genotype := fifelse(Genotype == "./.", NA_integer_,
                                 fifelse(Genotype == "0/0", 0L,
                                         fifelse(Genotype == "0/1", 1L,
                                                 fifelse(Genotype == "1/1", 2L, NA_integer_))))]
  # Add info score
  info <- fread(paste0("chr", chr, ".info.txt")) 
  geno_dat <- merge(geno_dat, info, by = "POS", all.x = TRUE)
  
  # Retain only  ID, Genotype and SNP
  geno_dat <- geno_dat %>%
    distinct(ID, SNP, Genotype)
  
  # Get the imputed data and reshape from wide to long
  imp_dat <- fread(paste0("ukb_chr_", chr, ".GP.txt")) %>%
    melt(id.vars = "POS", variable.name = "ID", value.name = "Genotype") %>%
    filter(Genotype != ".")
  imp_dat <- merge(imp_dat, info, by = "POS", all.x = TRUE)
  
  # Get map files
  imp_map <- imp_dat %>%
    distinct(SNP, POS)
  
  # Get the original genotypes (we can use this earlier to make pmissing (=1) larger for imputed)
  original_snps <- fread(paste0("/pub59/iasiimwe/ukbb_geno/ukb_snp_chr", chr, "_v2.bim"))$V2

  dat_all <- tibble()
  for (j in seq_along(imp_map$SNP)) {
    snp <- imp_map$SNP[[j]]
    geno_dat_snp <- geno_dat[SNP == snp]
    dat_snp <- imp_dat[SNP == snp]
    
    # Missingness percentage
    pmissing <- (length(dat_snp$ID) + geno_dat_snp %>% filter(is.na(Genotype)) %>% nrow()) / length(geno_dat_snp$ID)
    
    genotyped_status <- ifelse(snp %in% original_snps, "Genotyped", "Imputed")
    
    dat_all <- dat_all %>%
      bind_rows(tibble(SNP = snp, pmissing = pmissing, Genotyped = genotyped_status))
    
    message(paste0("Chr: ", chr, "; ", round(j * 100/length(imp_map$SNP), 3), "% complete!"))
  }
  write.table(dat_all, file = paste0("imputed_datasets/chr_", chr, "_pmissing.txt"), quote = FALSE, row.names = FALSE, sep = "\t")
}



# Perform traditional GWAS
# ------------------------
## Plink 2 defaults to analyzing all phenotypes/covariates 
# country age weight sex target_inr hiv_positive hiv_unknown sim_amio C1 C2 C3 C4 C5 C6 C7 C8 C9 C10
# Unless pheno_name or covar_name is used

# # Conduct QTL analysis
# for outcome in urate; do
# for nsample in 100 500 1000; do
# for i in 4; do
#  /pub59/iasiimwe/plink2 \
#   --bfile ukb_chr_${i} \
#   --pheno ukb_${nsample}_pheno.txt --pheno-name ${outcome} \
#   --covar ukb_${nsample}_covar.txt --covar-variance-standardize \
#   --glm \
#   --out ukb_${nsample}_chr_${i}
# done
# done
# done


# Perform GWAS with imputed data
# -----------------------------
# # Change to binary
# for chr in 4; do
# for i in {1..10}; do
# /pub59/iasiimwe/plink1.9/plink --file imputed_datasets/chr_${chr}_dat_${i} --make-bed --out imputed_datasets/chr_${chr}_dat_${i}
# done
# done
#
# for nsample in 100 500 1000; do
# for outcome in urate; do
# for chr in 4; do
# for i in {1..10}; do
#  /pub59/iasiimwe/plink2 \
#   --bfile imputed_datasets/chr_${chr}_dat_${i} \
#   --pheno ukb_${nsample}_pheno.txt --pheno-name ${outcome} \
#   --covar ukb_${nsample}_covar.txt --covar-variance-standardize \
#   --glm \
#   --out imputed_datasets/ukb_${nsample}_chr_${chr}_dat_${i}
# done
# done
# done
# done


# Get Manhattan data
library(tidyverse)
library(data.table)

# Rubin's rule function to pool estimates from multiple imputations
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


# Genotyped status
chr <- 4
geno_status <- fread(paste0("imputed_datasets/chr_", chr, "_pmissing.txt")) %>%
  select(SNP, Genotyped)
write_csv(geno_status, "chr4_geno_status.csv")

# Traditional GWAS
outcomes <- c("urate")
chrs <- c(4)
nsamples <- c(1000, 500, 100)
dat_all_nsample <- tibble()
for(nsample in nsamples) {
  dat_all_outcome <- tibble()
  for (outcome in outcomes) {
    dat_all <- tibble()
    for (chr in chrs) {
      info <- fread(paste0("chr", chr, ".info.txt")) %>%
        left_join(geno_status) %>%
        select(SNP, Genotyped)
      dat <- fread(paste0("ukb_", nsample, "_chr_", chr, ".", outcome, ".glm.linear")) %>%
        filter(TEST == "ADD") %>%
        select(CHR = "#CHROM", BP = POS, SNP = ID, BETA, SE, P) %>%
        mutate(Outcome = outcome, Nsample = nsample) %>%
        left_join(info)
      dat_all <- bind_rows(dat_all, dat)
    }
    dat_all_outcome <- bind_rows(dat_all_outcome, dat_all)
  }
  dat_all_nsample <- bind_rows(dat_all_nsample, dat_all_outcome)
}
write.table(dat_all_nsample, "Traditional_manhattan_ukb.txt", sep = " ", row.names = FALSE, quote = FALSE)

# Imputed
outcomes <- c("urate")
chrs <- c(4)
datasets <- c(1:10) # Also do for c(1:100)
nsamples <- c(1000, 500, 100)
dat_all_nsample <- tibble()
for(nsample in nsamples) {
  dat_all_outcome <- tibble()
  for (outcome in outcomes) {
    
    dat_all_chr <- tibble()
    for (chr in chrs) {
      info <- fread(paste0("chr", chr, ".info.txt")) %>%
        left_join(geno_status) %>%
        select(SNP, Genotyped)
      
      dat_all <- tibble()
      for (dataset in datasets) {
        dat <- fread(paste0("imputed_datasets/ukb_", nsample, "_chr_", chr, "_dat_", dataset, ".", outcome, ".glm.linear")) %>%
          filter(TEST == "ADD") %>%
          select(CHR = "#CHROM", BP = POS, SNP = ID, BETA, SE)
        colnames(dat) <- c("CHR", "BP", "SNP", paste0(c("BETA", "SE"), "_", dataset))
        if (dataset == 1) dat_all <- dat else dat_all <- left_join(dat_all, dat)
      }
      
      coef_estimates <- dat_all %>% select(contains("BETA"))
      coef_ses <- dat_all %>% select(contains("SE"))
      
      # Apply rubin_fn row-wise
      pooled_results <- pmap_dfr(
        list(split(coef_estimates, seq(nrow(coef_estimates))),
             split(coef_ses, seq(nrow(coef_ses)))),
        ~ rubin_fn(as.numeric(unlist(..1)), as.numeric(unlist(..2))) %>% as.list()
      )
      
      # Bind results to original dataset and save
      dat_all_pooled <- dat_all %>%
        select(SNP, CHR, BP) %>%
        bind_cols(pooled_results) %>%
        filter(!is.na(p)) %>%
        select(SNP, CHR, BP, BETA = pooled_est, SE = pooled_se, P = p_fmt) %>%
        mutate(Outcome = outcome, Nsample = nsample) %>%
        left_join(info)
      
      dat_all_chr <- bind_rows(dat_all_chr, dat_all_pooled)
    }
    dat_all_outcome <- bind_rows(dat_all_outcome, dat_all_chr)
  }
  dat_all_nsample <- bind_rows(dat_all_nsample, dat_all_outcome)
}
write.table(dat_all_nsample, paste0("Imputed_", length(datasets), "_datasets_manhattan_ukb.txt"), sep = " ", row.names = FALSE, quote = FALSE)



# Can do this locally 
# -------------------
library(tidyverse)
library(data.table)
library(qqman)
library(ggplotify) 
library(patchwork)
library(export)
library(ggrepel)
library(ggh4x)

# Read data
file_names <- c("Traditional_manhattan_ukb.txt",
                "Imputed_10_datasets_manhattan_ukb.txt",
                "Imputed_20_datasets_manhattan_ukb.txt",
                "Imputed_30_datasets_manhattan_ukb.txt") 
nice_names <- c("Single/Traditional", "Multiple (10 datasets)", 
                "Multiple (20 datasets)", "Multiple (30 datasets)")
dat <- tibble()
for (i in seq_along(file_names)) dat <- bind_rows(dat, mutate(fread(file_names[[i]]), Analysis = nice_names[[i]]))

# Get genotyped SNPs
geno_status <- read_csv("chr4_geno_status.csv", show_col_types = FALSE)

# Process
dat <- dat %>%
  filter(!is.na(P), Nsample != 100) %>%
  arrange(P) %>%
  select(-Genotyped) %>%
  left_join(geno_status) %>%
  mutate(Analysis = factor(Analysis, levels = c("Single/Traditional", "Multiple (10 datasets)", 
                                                "Multiple (20 datasets)", "Multiple (30 datasets)")),
         Nsample = factor(paste0("N = ", Nsample), levels = c("N = 500", "N = 1000")),
         # Nsample = factor(paste0("N = ", Nsample), levels = c("N = 100", "N = 500", "N = 1000")),
         Genotyped = if_else(is.na(Genotyped), "Imputed", Genotyped))

# Top SNPs (imputed and genotyped for N = 500)
top_snps <- dat %>%
  filter(Nsample == "N = 500") %>%
  arrange(P) %>%
  group_by(Genotyped, Analysis) %>%
  slice_head(n = 1) 
top_snps <- dat %>%
  filter(Nsample != "N = 100", SNP %in% top_snps$SNP) %>%
  mutate(
    label = paste0(SNP, "\n(β=", round(BETA, 1), "; P=", format(P, scientific = TRUE, digits = 2), ")"),
    # label = paste0(SNP, "\n(italic('\u03B2')=", round(BETA, 1),
    #                '*","~italic("P")=', format(P, scientific = TRUE, digits = 2), ")"), 
                   label = gsub("e-0", "E-", label))

# Compute cumulative BP position for plotting
dat <- dat %>%
  arrange(CHR, BP) %>%
  group_by(CHR) %>%
  mutate(index = row_number(),
         BP_spacing = index * 100) %>%  # Add space between SNPs in the same chromosome
  ungroup() %>%
  mutate(CHR_numeric = parse_number(as.character(CHR)),
         BP_cum = BP_spacing + (CHR_numeric - min(CHR_numeric)) * 5e4)  # Reduce gap between chromosomes
# top_n <- dat %>% filter(P < 5e-8) %>% group_by(Nsample) %>% count() # 500: 4, 1000: 243

# Create axis labels
axis_df <- dat %>%
  group_by(CHR) %>%
  summarize(center = mean(BP_cum))

#c("darkgray", "steelblue")

# Manhattan plot
pl <- ggplot(dat, aes(x = BP_cum, y = -log10(P))) +
  geom_point(alpha = 0.8, size = 1.2, colour = "darkgray") +
  geom_point(data = dat %>% filter(Genotyped == "Genotyped"), # black open circles around Genotyped SNPs
             shape = 21, fill = NA, color = "black", stroke = 0.5, size = 2) +
  geom_point(data = dat %>% filter(SNP %in% top_snps$SNP),
             colour = "red", alpha = 0.8, size = 1.2) +
  geom_text_repel(data = top_snps %>% left_join(dat), aes(label = label), parse = FALSE, size = 3, max.overlaps = Inf,
                  box.padding = 0.4, point.padding = 0.3, min.segment.length = 0, segment.color = "black",
                  segment.size = 0.2, nudge_y = 0.2) +
  facet_nested_wrap(vars(Analysis, Nsample), nrow = 2) +
  # scale_x_continuous(label = axis_df$CHR, breaks = axis_df$center) +
  geom_hline(yintercept = -log10(5e-8), linetype = "dashed", color = "red") +
  labs(x = "Chromosome 4 position", y = "-log10(P)", 
       title = "C. UK Biobank cohort") +
  theme_bw() +
  theme(legend.position = "none",
        strip.text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10, colour = "black"),
        # axis.text.x = element_text(face = "italic"),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0, vjust = 0, size = 12, face = "bold"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())
x <- 0.7
ggsave("Fig_3b_ukb.png", width = 20 * x, height = 10 * x)
# graph2ppt(pl, "Fig_3b_ukb", append = FALSE, width = 20 * x, height = 10 * x)







# Time spent on the real world data will need consideration because we tested a subset and only part of chromosome 4 - this shouldn't be a problem with parellization


# Covariate selection - use traditional GWAS
# -------------------
library(tidyverse)
library(data.table)
library(caret)
library(glmnet) # LASSO
library(randomForest) # RF

# Caret settings
fitControl <- trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 5)

# Dominant coding function
additive_coding_fn <- function(x) {
  x <- gsub(" ", "", x) 
  x[x == "00"] <- NA
  y <- table(x)
  y <- names(sort(-y)) # Sort in descending order, assume mutant type is the least abundant
  y <- paste0("^", y, "$")
  if (length(y) == 3) x <- gsub(y[[3]], "1", x) # Give mutant homozygotes "2"
  if (length(y) > 1) x <- gsub(y[[2]], "0.5", x) # Give heterozygotes "1" # So that we don't have to scale
  x <- gsub(y[[1]], "0", x) # Give homozygotes 0
  x[is.na(x)] <- "0" # Use mode imputation as dropping NAs leads to significant data loss
  return(as.numeric(x))
}

n <- 100
ml_methods <- c("glmnet", "rf")
outcomes <- c("urate")
chrs <- c(4)
nsamples <- c(1000, 500, 100)
dat_all_nsample <- tibble()
for(nsample in nsamples) {
  dat_all_outcome <- tibble()
  for (outcome in outcomes) {
    dat_all <- tibble()
    for (chr in chrs) {
      dat <- fread(paste0("ukb_", nsample, "_chr_", chr, ".", outcome, ".glm.linear")) %>%
        filter(TEST == "ADD") %>%
        arrange(P)
      dat[1:n, ] %>% 
        select(ID) %>%
        write.table(file = "topnsnps.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
      
      # Covariates
      covs <- tibble(fread(paste0("ukb_", nsample, "_covar.txt"))) %>%
        rename(ID = IID) %>%
        select(-FID)
      
      # Phenotype
      pheno <- tibble(fread(paste0("ukb_", nsample, "_pheno.txt")))[, c("IID", outcome)]
      colnames(pheno) <- c("ID", "y")
      
      # Sampled participants
      fread(paste0("ukb_", nsample, "_pheno.txt")) %>% 
        select(FID, IID) %>% 
        write.table(file = "temp_sampled.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
      
      # Get subset SNP data
      dat_path <- paste0("ukb_chr_", chr)
      temp_path <- "temp"
      system(paste0("/pub59/iasiimwe/plink1.9/plink --bfile ", dat_path, 
                    " --keep temp_sampled.txt --extract topnsnps.txt --recode tab --out ", temp_path))
      ped <- tibble(fread(paste0(temp_path, ".ped"), header = FALSE))
      map <- tibble(fread(paste0(temp_path, ".map"), header = FALSE))
      snps <- pull(map, V2)
      ped_start <- ped %>% select(V1:V6) # This selects FID, IID, father's ID, mother's ID, sex and phenotype data
      snp_columns <- colnames(ped)[!colnames(ped) %in% colnames(ped_start)]
      ped_snps <- ped %>%
        select(all_of(snp_columns)) 
      colnames(ped_snps) <- snps 
      ped_snps <- ped_snps %>%
        mutate_all(additive_coding_fn) %>%
        mutate(ID = ped_start$V2, .before = 1)

      # Put everything together
      ped_snps <- ped_snps %>%
        left_join(covs) %>%
        left_join(pheno) %>%
        select(-ID)
      
      # Perform the regression
      dat_method <- tibble()
      for (ml_method in ml_methods) {
        study <- ped_snps 
        set.seed(7)
        trained_mod <- train(y ~ ., data = study,
                             method = ml_method, 
                             trControl = fitControl, 
                             verbose = FALSE)
        
        # Variable importance
        top_covariates <- trained_mod %>% varImp()
        top_covariates <- tibble(Covariates = rownames(data.frame(top_covariates$importance)),
                                 Importance = top_covariates$importance$Overall) %>%
          arrange(desc(Importance)) %>%
          filter(str_detect(Covariates, "`|rs")) %>% # Remove non-SNP covariates
          filter(Importance >= 1) %>% # We will see how the correct or true SNPs were rated, if at all
          mutate(Outcome = outcome,
                 Chr = chr,
                 Method = ml_method,
                 Rank = row_number(),
                 Nsample = nsample)
        dat_method <- bind_rows(dat_method, top_covariates)
      }
      dat_all <- bind_rows(dat_all, dat_method)
    }
    dat_all_outcome <- bind_rows(dat_all_outcome, dat_all)
  }
  dat_all_nsample <- bind_rows(dat_all_nsample, dat_all_outcome)
}
write.table(dat_all_nsample, "ukb_ML_methods.txt", sep = " ", row.names = FALSE, quote = FALSE)


# Can do this locally 
# -------------------
library(tidyverse)
library(data.table)
library(qqman)
library(ggplotify) 
library(patchwork)
library(export)
library(ggrepel)
library(ggh4x)

# Read data and process
datML <- tibble(fread("ukb_ML_methods.txt")) %>%
  mutate(SNP = gsub("`", "", Covariates)) %>%
  filter(Rank < 6) %>% # Only the top 5
  select(-Covariates, -Importance, -Chr)
dat1 <- tibble(fread("Traditional_manhattan_ukb.txt")) %>%
  left_join(filter(datML, Method == "glmnet")) %>%
  filter(!is.na(P)) %>%
  mutate(Method = "glmnet")
dat2 <- tibble(fread("Traditional_manhattan_ukb.txt")) %>%
  left_join(filter(datML, Method == "rf")) %>%
  filter(!is.na(P)) %>%
  mutate(Method = "rf")

# Get genotyped SNPs
geno_status <- read_csv("chr4_geno_status.csv", show_col_types = FALSE)
# Process
dat <- dat1 %>%
  bind_rows(dat2) %>%
  filter(Nsample != 100) %>%
  arrange(P) %>%
  select(-Genotyped) %>%
  left_join(geno_status) %>%
  mutate(Analysis = case_when(Method == "rf" ~ "Random Forest", Method == "glmnet" ~ "Penalized regression"),
         Analysis = factor(Analysis, levels = c("Random Forest", "Penalized regression")),
         Nsample = factor(paste0("N = ", Nsample), levels = c("N = 500", "N = 1000")),
         Genotyped = if_else(is.na(Genotyped), "Imputed", Genotyped))

# Top SNPs
top_snps <- dat %>%
  filter(!is.na(Rank)) %>%
  mutate(label = paste0(SNP, " (", Rank, ")"))

# Compute cumulative BP position for plotting
dat <- dat %>%
  arrange(CHR, BP) %>%
  group_by(CHR) %>%
  mutate(index = row_number(),
         BP_spacing = index * 100) %>%  # Add space between SNPs in the same chromosome
  ungroup() %>%
  mutate(CHR_numeric = parse_number(as.character(CHR)),
         BP_cum = BP_spacing + (CHR_numeric - min(CHR_numeric)) * 5e4)  # Reduce gap between chromosomes
# top_n <- dat %>% filter(P < 5e-8) %>% group_by(Nsample) %>% count() # 500: 4, 1000: 243

# Create axis labels
axis_df <- dat %>%
  group_by(CHR) %>%
  summarize(center = mean(BP_cum))

#c("darkgray", "steelblue")

# Manhattan plot
pl <- ggplot(dat, aes(x = BP_cum, y = -log10(P))) +
  geom_point(alpha = 0.8, size = 1.2, colour = "darkgray") +
  geom_point(data = dat %>% filter(Genotyped == "Genotyped"), # black open circles around Genotyped SNPs
             shape = 21, fill = NA, color = "black", stroke = 0.5, size = 2) +
  geom_point(data = dat %>% filter(SNP %in% top_snps$SNP, !is.na(Rank)),
             colour = "red", alpha = 0.8, size = 1.2) +
  geom_text_repel(data = top_snps %>% left_join(dat), aes(label = label), parse = FALSE, size = 3, max.overlaps = Inf,
                  box.padding = 0.4, point.padding = 0.3, min.segment.length = 0, segment.color = "black",
                  segment.size = 0.2, nudge_y = 0.2) +
  facet_nested_wrap(vars(Analysis, Nsample), nrow = 2) +
  # scale_x_continuous(label = axis_df$CHR, breaks = axis_df$center) +
  geom_hline(yintercept = -log10(5e-8), linetype = "dashed", color = "red") +
  labs(x = "Chromosome 4 position", y = "-log10(P)", 
       title = "C. UK Biobank cohort") +
  theme_bw() +
  theme(legend.position = "none",
        strip.text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10, colour = "black"),
        # axis.text.x = element_text(face = "italic"),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0, vjust = 0, size = 12, face = "bold"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())
x <- 0.7
ggsave("Fig_6b_ukb.png", width = 10 * x, height = 10 * x)
# graph2ppt(pl, "Fig_6b_ukb", append = FALSE, width = 10 * x, height = 10 * x)
