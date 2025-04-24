# Start with per-sample QC'd data
# -------------------------------
# cd /pub59/iasiimwe/TB/warpath; source /pub59/iasiimwe/miniconda3/bin/activate base

# zcat /pub59/iasiimwe/Tractor/Imputed/Black_original/chr10_qc_bgzip.vcf.gz | head -18
##fileformat=VCFv4.1
##filedate=2021.8.25
##contig=<ID=10>
##pipeline=michigan-imputationserver-1.5.7
##imputation=minimac4-1.0.2
##phasing=eagle-2.4
##r2Filter=0.0
##INFO=<ID=AF,Number=1,Type=Float,Description="Estimated Alternate Allele Frequency">
##INFO=<ID=MAF,Number=1,Type=Float,Description="Estimated Minor Allele Frequency">
##INFO=<ID=R2,Number=1,Type=Float,Description="Estimated Imputation Accuracy (R-square)">
##INFO=<ID=ER2,Number=1,Type=Float,Description="Empirical (Leave-One-Out) R-square (available only for genotyped variants)">
##INFO=<ID=IMPUTED,Number=0,Type=Flag,Description="Marker was imputed but NOT genotyped">
##INFO=<ID=TYPED,Number=0,Type=Flag,Description="Marker was genotyped AND imputed">
##INFO=<ID=TYPED_ONLY,Number=0,Type=Flag,Description="Marker was genotyped but NOT imputed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DS,Number=1,Type=Float,Description="Estimated Alternate Allele Dosage : [P(0/1)+2*P(1/1)]">
##FORMAT=<ID=HDS,Number=2,Type=Float,Description="Estimated Haploid Alternate Allele Dosage ">
##FORMAT=<ID=GP,Number=3,Type=Float,Description="Estimated Posterior Probabilities for Genotypes 0/0, 0/1 and 1/1 "

# # Get Plink files (standard analysis)
  # Base-pair positions for CYP2C9: 96698415..96749848 (Assembly: GRCh37.p13 (GCF_000001405.25), https://www.ncbi.nlm.nih.gov/gene/1559)
  # BP for VKORC1: 31102175..31106118 (Assembly: GRCh37.p13 (GCF_000001405.25), https://www.ncbi.nlm.nih.gov/gene/79001)
# Add 10,000 BPs on either side i.e. 
   # CYP2C9: 96688415..96759848
   # VKORC1: 31092175..31116118

# for i in 10 16; do
# if [ "$i" -eq 10 ]; then bp1=96688415; bp2=96759848; else bp1=31092175; bp2=31116118; fi
# /pub59/iasiimwe/plink2 \
#  --vcf /pub59/iasiimwe/Tractor/Imputed/Black_original/chr${i}_qc_bgzip.vcf.gz \
#  --chr ${i} --from-bp ${bp1} --to-bp ${bp2} \
#  --geno 0.05 --maf 0.01 --hwe 0.000001 --make-bed \
#  --out warpath_chr_${i}
# done

# # Genotyped
# for i in 10 16; do
# if [ "$i" -eq 10 ]; then bp1=96688415; bp2=96759848; else bp1=31092175; bp2=31116118; fi
# /pub59/iasiimwe/plink1.9/plink \
#  --bfile /pub59/iasiimwe/WARPATH/WARPATH \
#  --chr ${i} --from-bp ${bp1} --to-bp ${bp2} \
#  --geno 0.05 --maf 0.01 --hwe 0.000001 \
#  --keep samples_to_keep.txt \
#  --recode tab --out original_warpath_chr_${i}
# done

# # Extract genotypes (best guess, GT) and probabilities (GP)
# for i in 10 16; do
# # Get genotypes
# /pub59/iasiimwe/vcftools_0.1.13/bin/vcftools \
# --gzvcf /pub59/iasiimwe/Tractor/Imputed/Black_original/chr${i}_qc_bgzip.vcf.gz \
# --out warpath_chr_${i} \
# --extract-FORMAT-info GT
# 
# # Get probabilities
# /pub59/iasiimwe/vcftools_0.1.13/bin/vcftools \
#   --gzvcf /pub59/iasiimwe/Tractor/Imputed/Black_original/chr${i}_qc_bgzip.vcf.gz \
#   --out warpath_chr_${i} \
#   --extract-FORMAT-info GP
# done

# # For GP, subset to only the CYP2C9 and VKORC1 BPs
# -------------------------------------------------
library(tidyverse)
library(data.table)

for (i in c(10, 16)) {
  # Read bim file to get the applicable SNPs
  bim <- fread(paste0("warpath_chr_", i, ".bim")) %>%
    select(SNP = V2, POS = V4)
  
  # Get info scores
  info <- fread(paste0("/pub59/iasiimwe/INFO_files/Discovery/chr", i, ".info.gz")) %>%
    filter(SNP %in% bim$SNP) %>% # MAF >= 0.01, 
    select(SNP,  MAF, AvgCall, Rsq, Genotyped) %>%
    rowwise() %>%
    mutate(POS = strsplit(SNP, ":")[[1]][2]) %>%
    arrange(POS, desc(Rsq)) %>%
    group_by(POS) %>%
    slice_head(n = 1) 
  fwrite(info, paste0("chr", i, ".info.txt"))
  
  # Genotype probabilities
  imp_dat <- fread(paste0("warpath_chr_", i, ".GP.FORMAT")) %>%
    select(-CHROM) %>%
    filter(POS %in% bim$POS)
  fwrite(imp_dat, paste0("warpath_chr_", i, ".GP.txt"))
}


# Write the imputed datasets
# --------------------------
library(tidyverse)
library(data.table)

# Add folders if they don't exist
if (!dir.exists("imputed_datasets")) dir.create("imputed_datasets")

# Obtain original data

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

# Replace NA's function
replace_non_na <- function(x, y) {
  # Identify the two most common genotype strings in `y`
  most_common <- names(sort(table(y), decreasing = TRUE))[1:2]
  
  # Extract unique alleles from the most common genotypes
  alleleA <- unique(unlist(str_split(most_common[1], " "))) # maps to A
  alleleB <- unique(unlist(str_split(most_common[2], " ")))
  alleleB <- setdiff(alleleB, alleleA) # maps to T
  
  if (length(alleleB) == 0) { # Heterozygotes more than homozygotes
    alleleA <- unique(unlist(str_split(most_common[2], " "))) # maps to A
    alleleB <- unique(unlist(str_split(most_common[1], " ")))
    alleleB <- setdiff(alleleB, alleleA) # maps to T
  }
  
  # Defensive check
  if (length(alleleA) != 1 || length(alleleB) != 1) {
    stop("Ambiguous or missing allele mappings. Check input data.")
  }
  
  # Recode alleles
  y_recoded <- y %>%
    str_replace_all(alleleA, "1") %>%
    str_replace_all(alleleB, "2") %>%
    str_replace_all("1", "A") %>%
    str_replace_all("2", "T")
  
  # Enforce consistent ordering
  y_recoded[y_recoded == "T A"] <- "A T"
  
  # Replace NA values in x with recoded y
  x_filled <- ifelse(is.na(x), y_recoded, x)
  
  return(x_filled)
}

# Do per chromosome
chrs <- c(10, 16)
for (chr in chrs) {
  # Get and process the original data
  ped <- tibble(fread(paste0("original_warpath_chr_", chr, ".ped"), header = FALSE))
  map <- tibble(fread(paste0("original_warpath_chr_", chr, ".map"), header = FALSE))
  snps <- pull(map, V2)
  ped_start <- ped %>% select(V1:V6) # This selects FID, IID, father's ID, mother's ID, sex and phenotype data
  snp_columns <- colnames(ped)[!colnames(ped) %in% colnames(ped_start)]
  ped_snps <- ped %>%
    select(all_of(snp_columns)) 
  colnames(ped_snps) <- snps 
  ped_snps <- ped_snps %>%
    mutate(ID = paste0(ped_start$V1, "_", ped_start$V2), .before = 0)
  map <- map %>% select(SNP = V2, POS = V4) %>% group_by(POS) %>% slice_head(n = 1)

  # Get the imputed data
  imp_dat <- fread(paste0("warpath_chr_", chr, ".GP.txt")) 

  # Reshape from wide to long (gather samples under an "ID" column)
  imp_dat <- melt(imp_dat, id.vars = "POS", variable.name = "ID", value.name = "Genotype") %>%
    filter(Genotype != ".")
  
  # Add info score
  info <- fread(paste0("chr", chr, ".info.txt")) 
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
    select(ID, Genotype, SNP, R2 = Rsq, MAF, Genotyped, POS)
  
  nsim <- 100
  set.seed(7)
  for (i in 1:nsim) {
    
    dat_all <- ped_start
    for (j in seq_along(imp_map$SNP)) {
      snp <- imp_map$SNP[[j]]
      dat_snp <- imp_dat[SNP == snp]
      
      # Missingness percentage
      if (unique(dat_snp$Genotyped) == "Imputed") {
        pmissing <- 1
      } else {
        col <- map %>% filter(POS == unique(dat_snp$POS)) %>% pull(SNP)
        if (length(col) == 0) pmissing <- 1 else {
          geno_dat_snp <- ped_snps %>% select(ID, contains(col))
          colnames(geno_dat_snp) <- c("ID", "osnp")
          pmissing <- (geno_dat_snp %>% filter(is.na(osnp)) %>% nrow())/nrow(geno_dat_snp)
        }
      }

      # Add coded genotypes
      dat_snp[, Coded := additive_coding_fn(Genotype, R2, pmissing)]
      set(dat_snp, which(is.na(dat_snp$Coded)), "Coded", "0 0")
      dat_snp <- dat_snp[, .SD[1], by = ID]
      
      if (unique(dat_snp$Genotyped) == "Genotyped") {
        col <- map %>% filter(POS == unique(dat_snp$POS)) %>% pull(SNP)
        if (length(col) == 1) {
          geno_dat_snp <- ped_snps %>% select(ID, contains(col))
          colnames(geno_dat_snp) <- c("ID", "osnp")
          dat_snp$Coded <- replace_non_na(dat_snp$Coded, geno_dat_snp$osnp)
        }
        }
      
      # Add to dat_all
      dat_all <- merge(dat_all, dat_snp[, .(ID, Coded)], by = "ID", all.x = TRUE) # Imputed data
      setnames(dat_all, "Coded", snp)
      
      message(paste0("Chr: ", chr, "; dataset ", i, " of ", nsim, "; ", round(j * 100/length(imp_map$SNP), 3), "% complete!"))
    }
    
    # Save map and ped files
    
    write.table(imp_map, file = paste0("imputed_datasets/chr_", chr, "_dat_", i, ".map"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
    write.table(dat_all, file = paste0("imputed_datasets/chr_", chr, "_dat_", i, ".ped"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  }
}


# Start with dose and then other outcomes e.g. s-warfarin/R-warfarin

# Get covariates and phenotypes (R) # https://pmc.ncbi.nlm.nih.gov/articles/PMC9537548/
# ---------------------------------
library(tidyverse)
library(data.table)

dat <- fread("/pub59/iasiimwe/Tractor/Imputed/PK/WARPATH_PK_covariates.sample") 
dat <- dat[-1, ] # Remove first row

# Covariates
covs <- dat %>%
  mutate(FID = 0, IID = paste0(ID_1, "_", ID_2)) %>%
  select(FID, IID, age = Age_years, sex = GenderMale, weight = Weight_kg, target_inr = Target_INRThree, 
         efv = Efavirenz_EFVYes, sim_amio = Simvastatin_amiodaroneYes, C1:C2)
write.table(covs, "warpath_covariates.txt", sep = " ", row.names = FALSE, quote = FALSE)

# Covariates imputed
covs <- dat %>%
  mutate(FID = paste0(ID_1, "_", ID_2), IID = FID) %>%
  select(FID, IID, age = Age_years, sex = GenderMale, weight = Weight_kg, target_inr = Target_INRThree, 
         efv = Efavirenz_EFVYes, sim_amio = Simvastatin_amiodaroneYes, C1:C2)
write.table(covs, "warpath_covariates_imputed.txt", sep = " ", row.names = FALSE, quote = FALSE)

# Pheno
pheno <- dat %>%
  mutate(FID = 0, IID = paste0(ID_1, "_", ID_2)) %>%
  select(FID, IID, Ln_dose_mg_week:RS_all_met_ratio)
write.table(pheno, "warpath_phenotypes.txt", sep = " ", row.names = FALSE, quote = FALSE)

# Pheno imputed
pheno <- dat %>%
  mutate(FID = paste0(ID_1, "_", ID_2), IID = FID) %>%
  select(FID, IID, Ln_dose_mg_week:RS_all_met_ratio)
write.table(pheno, "warpath_phenotypes_imputed.txt", sep = " ", row.names = FALSE, quote = FALSE)


# Perform traditional GWAS
# ------------------------
## Plink 2 defaults to analyzing all phenotypes/covariates 
  # country age weight sex target_inr hiv_positive hiv_unknown sim_amio C1 C2 C3 C4 C5 C6 C7 C8 C9 C10
  # Unless pheno_name or covar_name is used

# # Conduct QTL analysis
# for outcome in Ln_dose_mg_week S_R_ratio; do
# for i in 10 16; do
#  /pub59/iasiimwe/plink2 \
#   --bfile warpath_chr_${i} \
#   --pheno warpath_phenotypes.txt --pheno-name ${outcome} \
#   --covar warpath_covariates.txt --covar-variance-standardize \
#   --glm \
#   --out warpath_chr_${i}
# done
# done

# Perform GWAS with imputed data
# -----------------------------
# # Change to binary
# for chr in 10 16; do
# for i in {1..100}; do
# /pub59/iasiimwe/plink1.9/plink --file imputed_datasets/chr_${chr}_dat_${i} --make-bed --out imputed_datasets/chr_${chr}_dat_${i}
# done
# done
#
# for outcome in Ln_dose_mg_week S_R_ratio; do
# for chr in 10 16; do
# for i in {1..100}; do
#  /pub59/iasiimwe/plink2 \
#   --bfile imputed_datasets/chr_${chr}_dat_${i} \
#   --pheno warpath_phenotypes_imputed.txt --pheno-name ${outcome} \
#   --covar warpath_covariates_imputed.txt --covar-variance-standardize \
#   --glm \
#   --out imputed_datasets/warpath_chr_${chr}_dat_${i}
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

# Traditional GWAS
outcomes <- c("Ln_dose_mg_week", "S_R_ratio")
chrs <- c(10, 16)
dat_all_outcome <- tibble()
for (outcome in outcomes) {
  dat_all <- tibble()
  for (chr in chrs) {
   info <- fread(paste0("chr", chr, ".info.txt")) %>%
      select(BP = POS, Genotyped)
   dat <- fread(paste0("warpath_chr_", chr, ".", outcome, ".glm.linear")) %>%
     filter(TEST == "ADD") %>%
     select(CHR = "#CHROM", BP = POS, SNP = ID, BETA, SE, P) %>%
     mutate(Outcome = outcome) %>%
     left_join(info)
   dat_all <- bind_rows(dat_all, dat)
  }
  dat_all_outcome <- bind_rows(dat_all_outcome, dat_all)
}
write.table(dat_all_outcome, "Traditional_manhattan.txt", sep = " ", row.names = FALSE, quote = FALSE)

# Imputed
datasets <- c(1:10) # Also do for c(1:100)
dat_all_outcome <- tibble()
for (outcome in outcomes) {
  
  dat_all_chr <- tibble()
  for (chr in chrs) {
    info <- fread(paste0("chr", chr, ".info.txt")) %>%
      select(BP = POS, Genotyped)

    dat_all <- tibble()
    for (dataset in datasets) {
      dat <- fread(paste0("imputed_datasets/warpath_chr_", chr, "_dat_", dataset, ".", outcome, ".glm.linear")) %>%
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
      mutate(Outcome = outcome) %>%
      left_join(info)
    
    dat_all_chr <- bind_rows(dat_all_chr, dat_all_pooled)
  }
  dat_all_outcome <- bind_rows(dat_all_outcome, dat_all_chr)
}
write.table(dat_all_outcome, paste0("Imputed_", length(datasets), "_datasets_manhattan.txt"), sep = " ", row.names = FALSE, quote = FALSE)


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

known_snps <- tribble(
  ~SNP,                        ~rsID,                     ~Common_name,         ~CHR, ~BP,
  "10:96702047:C:T",          "rs1799853",              "CYP2C9*2",           10,    96702047,
  "10:96741053:A:C",          "rs1057910",              "CYP2C9*3",           10,    96741053,
  "10:96741058:C:G",          "rs28371686",             "CYP2C9*5",           10,    96741058,
  "10:96709038:GA:G",         "rs9332131 (GA/G)",       "CYP2C9*6",           10,    96709038,
  "10:96702066:G:A",          "rs7900194 (G/A)",        "CYP2C9*8",           10,    96702066,
  "10:96708974:A:G",          "rs2256871 (A/G)",        "CYP2C9*9",           10,    96708974,
  "10:96740981:C:T",          "rs28371685 (C/T)",       "CYP2C9*11",          10,    96740981,
  "10:96405502:G:A",          "rs12777823 (G/A)",       NA,                   10,    96405502,
  # "16:31102321:C:T",          "rs7294 (G/A)",           "VKORC1 3730G>A",     16,    31102321,
  # "16:31103796:A:G",          "rs2359612 (C/T)",        "VKORC1 2255C>T",     16,    31103796,
  # "16:31104509:C:G",          "rs8050894 (G/C)",        "VKORC1 1542G>C",     16,    31104509,
  # "16:31104878:G:A",          "rs9934438 (C/T)",        "VKORC1 1173C>T",     16,    31104878,
  # "16:31105554:A:C",          "rs2884737 (T/G)",        "VKORC1 497T>G",      16,    31105554,
  "16:31107689:C:T",          "rs9923231 (G/A)",        "VKORC1 -1639G>A",    16,    31107689 #,
  #"19:15990431:C:T",          "rs2108622 (C/T)",        "CYP4F2*3",           19,    15990431
)

# Read data
file_names <- c("Traditional_manhattan.txt",
                "Imputed_10_datasets_manhattan.txt") # Focus on 10 as no additional benefit with 100
# nice_names <- c("Single/Traditional", "Multiple (10 datasets)", "Multiple (100 datasets)")
nice_names <- c("Single Imputation (Traditional approach)", "Multiple Imputation")
dat <- tibble()
for (i in seq_along(file_names)) dat <- bind_rows(dat, mutate(fread(file_names[[i]]), Analysis = nice_names[[i]]))

# Process
dat <- dat %>%
  left_join(select(known_snps, SNP, Common_name)) %>%
  filter(!is.na(P)) %>%
  arrange(P) %>%
  mutate(
    # CHR = ifelse(CHR == 10, 
    #                   paste0(CHR, " (CYP2C9)"),
    #                   paste0(CHR, " (VKORC1)")),
         CHR = factor(CHR, levels = sort(unique(CHR))),
         Outcome = case_when(Outcome == "Ln_dose_mg_week" ~ "Weekly warfarin dose",
                             Outcome == "S_R_ratio" ~ "S:R warfarin concentration ratios"),
         Analysis = factor(Analysis, levels = c("Single Imputation (Traditional approach)", "Multiple Imputation")),
         Common_name = if_else(is.na(Common_name) | P > 1e-2, "", Common_name),
         Common_name = if_else(BP == 96758308, "rs58800757", Common_name),
         Common_name_label = if_else(Common_name != "", 
                                     paste0("italic('", Common_name, "')"), 
                                     NA_character_),
         SNP_group = if_else(Common_name != "", "highlight", as.character(CHR)))
# dat %>% filter(Common_name != "") 

# Compute cumulative BP position for plotting
dat <- dat %>%
  arrange(CHR, BP) %>%
  group_by(CHR) %>%
  mutate(index = row_number(),
         BP_spacing = index * 100) %>%  # Add space between SNPs in the same chromosome
  ungroup() %>%
  mutate(CHR_numeric = parse_number(as.character(CHR)),
         BP_cum = BP_spacing + (CHR_numeric - min(CHR_numeric)) * 5e4)  # Reduce gap between chromosomes
# top_n <- dat %>% filter(P < 5e-8) %>% group_by(Outcome) %>% count()

# Create axis labels
axis_df <- dat %>%
  group_by(CHR) %>%
  summarize(center = mean(BP_cum))

# Colour pallette
chr_levels <- unique(as.character(dat$CHR))
palette_chr <- setNames(rep(c("darkgray", "steelblue"), length.out = length(chr_levels)), chr_levels)
custom_colors <- c(palette_chr, "highlight" = "red")

# Manhattan plot
pl <- ggplot(dat, aes(x = BP_cum, y = -log10(P))) +
  geom_point(aes(color = SNP_group), alpha = 0.8, size = 1.2) +
  geom_point(data = dat %>% filter(Genotyped == "Genotyped"), # black open circles around Genotyped SNPs
             shape = 21, fill = NA, color = "black", stroke = 0.5, size = 2) +
  geom_point(data = dat %>% filter(Common_name != ""),
             colour = "red", alpha = 0.8, size = 1.2) +
  scale_color_manual(values = custom_colors) +
  geom_text_repel(data = dat %>% filter(Common_name != ""),
                  aes(label = Common_name_label),
                  parse = TRUE,
                  size = 3,
                  max.overlaps = Inf,
                  box.padding = 0.4,
                  point.padding = 0.3,
                  min.segment.length = 0,
                  segment.color = "black",
                  segment.size = 0.2,
                  nudge_y = 0.2) +
  facet_nested_wrap(vars(Analysis, Outcome), nrow = 3) +
  scale_x_continuous(label = axis_df$CHR, breaks = axis_df$center) +
  geom_hline(yintercept = -log10(5e-8), linetype = "dashed", color = "red") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  labs(x = "Chromosome", y = "-log10(P)", 
       title = "B. War-PATH cohort") +
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
graph2ppt(pl, "Fig_3a", append = FALSE, width = 10 * x, height = 10 * x)

# Able to detect known associations (including imputed)
# Fewer false positives
  # rs58800757 in LD with *8, imputed and gave a false association - also many false signals around VKORC1
  # Message is that multiple imputation is more likely to detect a true/causal SNP (*8 was also imputed)

dat %>%
  filter(P < 5e-8) %>%
  group_by(CHR, Outcome, Analysis) %>%
  count()


# Check in the Manhattan plots their p-values
   # dat %>% filter(Common_name != "")
# Check in info scores if they were imputed 
# If imputed, variability should be captured (at the moment, it isn't as p-values are similar)
# Check in dat and see if the p-values are the same
# Can we include those with lower MAF because we have uncertainty accounted for?



known_snps <- known_snps %>% select(SNP, Common_name)
fread("chr10.info.txt") %>% left_join(known_snps) %>% filter(!is.na(Common_name))
# Joining with `by = join_by(SNP)`
# SNP     Rsq Genotyped      POS Common_name
# 1:  10:96702066:G:A 0.99910   Imputed 96702066    CYP2C9*8
# 2:  10:96708974:A:G 0.99696 Genotyped 96708974    CYP2C9*9
# 3: 10:96709038:GA:G 0.99542   Imputed 96709038    CYP2C9*6
# 4:  10:96740981:C:T 0.96852 Genotyped 96740981   CYP2C9*11
fread("chr16.info.txt") %>% left_join(known_snps) %>% filter(!is.na(Common_name))
# SNP     Rsq Genotyped      POS     Common_name
# 1: 16:31107689:C:T 0.99969 Genotyped 31107689 VKORC1 -1639G>A



# # LD between CYP2C9*8 and rs58800757
# /pub59/iasiimwe/plink1.9/plink --bfile warpath_chr_10 --ld 10:96758308:G:A 10:96702066:G:A --r2 --out ld_rs58800757_rs7900194_8 # R-sq = 0.118372
# # LD between CYP2C9*9 and rs58800757
# /pub59/iasiimwe/plink1.9/plink --bfile warpath_chr_10 --ld 10:96758308:G:A 10:96708974:A:G --r2 --out ld_rs58800757_rs2256871_9 # R-sq = 0.0261663
# # LD between CYP2C9*11 and rs58800757
# /pub59/iasiimwe/plink1.9/plink --bfile warpath_chr_10 --ld 10:96758308:G:A 10:96740981:C:T --r2 --out ld_rs58800757_rs28371685_11 # R-sq = 0.0321845



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
  x[x == "0 0"] <- NA
  y <- table(x)
  y <- names(sort(-y)) # Sort in descending order, assume mutant type is the least abundant
  y <- paste0("^", y, "$")
  if (length(y) == 3) x <- gsub(y[[3]], "1", x) # Give mutant homozygotes "2"
  if (length(y) > 1) x <- gsub(y[[2]], "0.5", x) # Give heterozygotes "1" # So that we don't have to scale
  x <- gsub(y[[1]], "0", x) # Give homozygotes 0
  return(as.numeric(x))
}

n <- 100
ml_methods <- c("glmnet", "rf")
# Traditional GWAS
outcomes <- c("Ln_dose_mg_week", "S_R_ratio")
chrs <- c(10, 16)
dat_all_outcome <- tibble()
for (outcome in outcomes) {
  dat_all <- tibble()
  for (chr in chrs) {
    dat <- fread(paste0("warpath_chr_", chr, ".", outcome, ".glm.linear")) %>%
      filter(TEST == "ADD") %>%
      arrange(P)
    dat[1:n, ] %>% 
      select(ID) %>%
      write.table(file = "topnsnps.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    # Get subset SNP data
    dat_path <- paste0("warpath_chr_", chr)
    temp_path <- "temp"
    system(paste0("/pub59/iasiimwe/plink1.9/plink --bfile ", dat_path, 
                  " --extract topnsnps.txt --recode tab --out ", temp_path))
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

    # Covariates
    covs <- fread("warpath_covariates.txt") %>%
      rename(ID = IID) %>%
      select(-FID)
    
    # Phenotype
    pheno <- tibble(fread("warpath_phenotypes.txt"))[, c("IID", outcome)]
    colnames(pheno) <- c("ID", "y")
    
    # Put everything together
    ped_snps <- ped_snps %>%
      left_join(covs) %>%
      left_join(pheno) %>%
      select(-ID)
    
    # Perform the regression
    dat_method <- tibble()
    for (ml_method in ml_methods) {
      study <- ped_snps  %>%
        filter(!is.na(y))
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
        filter(str_detect(Covariates, "`")) %>% # Remove non-SNP covariates
        filter(Importance >= 1) %>% # We will see how the correct or true SNPs were rated, if at all
        mutate(Outcome = outcome,
               Chr = chr,
               Method = ml_method,
               Rank = row_number())
      dat_method <- bind_rows(dat_method, top_covariates)
    }
    dat_all <- bind_rows(dat_all, dat_method)
  }
  dat_all_outcome <- bind_rows(dat_all_outcome, dat_all)
}
write.table(dat_all_outcome, "WarPATH_ML_methods.txt", sep = " ", row.names = FALSE, quote = FALSE)




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

known_snps <- tribble(
  ~SNP,                        ~rsID,                     ~Common_name,         ~CHR, ~BP,
  "10:96702047:C:T",          "rs1799853",              "CYP2C9*2",           10,    96702047,
  "10:96741053:A:C",          "rs1057910",              "CYP2C9*3",           10,    96741053,
  "10:96741058:C:G",          "rs28371686",             "CYP2C9*5",           10,    96741058,
  "10:96709038:GA:G",         "rs9332131 (GA/G)",       "CYP2C9*6",           10,    96709038,
  "10:96702066:G:A",          "rs7900194 (G/A)",        "CYP2C9*8",           10,    96702066,
  "10:96708974:A:G",          "rs2256871 (A/G)",        "CYP2C9*9",           10,    96708974,
  "10:96740981:C:T",          "rs28371685 (C/T)",       "CYP2C9*11",          10,    96740981,
  "10:96405502:G:A",          "rs12777823 (G/A)",       NA,                   10,    96405502,
  "16:31107689:C:T",          "rs9923231 (G/A)",        "VKORC1 -1639G>A",    16,    31107689 
)

# Read data and process
datML <- tibble(fread("WarPATH_ML_methods.txt")) %>%
  mutate(SNP = gsub("`", "", Covariates)) %>%
  filter(Rank < 6) %>% # Only the top 5
  select(-Covariates, -Importance, -Chr)

dat1 <- tibble(fread("Traditional_manhattan.txt")) %>%
  left_join(filter(datML, Method == "glmnet")) %>%
  left_join(select(known_snps, SNP, Common_name)) %>%
  filter(!is.na(P)) %>%
  mutate(Method = "glmnet")
dat2 <- tibble(fread("Traditional_manhattan.txt")) %>%
  left_join(filter(datML, Method == "rf")) %>%
  left_join(select(known_snps, SNP, Common_name)) %>%
  filter(!is.na(P)) %>%
  mutate(Method = "rf")
dat <- dat1 %>%
  bind_rows(dat2) %>%
  arrange(P) %>%
  mutate(
    CHR = factor(CHR, levels = sort(unique(CHR))),
    Outcome = case_when(Outcome == "Ln_dose_mg_week" ~ "Weekly warfarin dose",
                        Outcome == "S_R_ratio" ~ "S:R warfarin concentration ratios"),
    Analysis = case_when(Method == "rf" ~ "Random Forest", Method == "glmnet" ~ "Penalized regression"),
    Analysis = factor(Analysis, levels = c("Random Forest", "Penalized regression")),
    Common_name = if_else(is.na(Common_name) | P > 1e-2, "", Common_name),
    Common_name = if_else(is.na(Rank), Common_name, paste0(Common_name, " (", Rank, ")")),
    Common_name_label = if_else(Common_name != "", 
                                paste0("italic('", Common_name, "')"), 
                                NA_character_),
    SNP_group = if_else(Common_name != "", "highlight", as.character(CHR)))
# dat %>% filter(Common_name != "") 

# Compute cumulative BP position for plotting
dat <- dat %>%
  arrange(CHR, BP) %>%
  group_by(CHR) %>%
  mutate(index = row_number(),
         BP_spacing = index * 100) %>%  # Add space between SNPs in the same chromosome
  ungroup() %>%
  mutate(CHR_numeric = parse_number(as.character(CHR)),
         BP_cum = BP_spacing + (CHR_numeric - min(CHR_numeric)) * 5e4)  # Reduce gap between chromosomes
# top_n <- dat %>% filter(P < 5e-8) %>% group_by(Outcome) %>% count()

# Create axis labels
axis_df <- dat %>%
  group_by(CHR) %>%
  summarize(center = mean(BP_cum))

# Colour pallette
chr_levels <- unique(as.character(dat$CHR))
palette_chr <- setNames(rep(c("darkgray", "steelblue"), length.out = length(chr_levels)), chr_levels)
custom_colors <- c(palette_chr, "highlight" = "red")

# Manhattan plot
pl <- ggplot(dat, aes(x = BP_cum, y = -log10(P))) +
  geom_point(aes(color = SNP_group), alpha = 0.8, size = 1.2) +
  geom_point(data = dat %>% filter(Genotyped == "Genotyped"), # black open circles around Genotyped SNPs
             shape = 21, fill = NA, color = "black", stroke = 0.5, size = 2) +
  geom_point(data = dat %>% filter(Common_name != ""),
             colour = "red", alpha = 0.8, size = 1.2) +
  scale_color_manual(values = custom_colors) +
  geom_text_repel(data = dat %>% filter(Common_name != ""),
                  aes(label = Common_name_label),
                  parse = TRUE,
                  size = 3,
                  max.overlaps = Inf,
                  box.padding = 0.4,
                  point.padding = 0.3,
                  min.segment.length = 0,
                  segment.color = "black",
                  segment.size = 0.2,
                  nudge_y = 0.2) +
  facet_nested_wrap(vars(Analysis, Outcome), nrow = 3) +
  scale_x_continuous(label = axis_df$CHR, breaks = axis_df$center) +
  geom_hline(yintercept = -log10(5e-8), linetype = "dashed", color = "red") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  labs(x = "Chromosome", y = "-log10(P)", 
       title = "B. War-PATH cohort") +
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
graph2ppt(pl, "Fig_6a", append = FALSE, width = 10 * x, height = 10 * x)





