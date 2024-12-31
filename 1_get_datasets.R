# Load libraries
library(data.table)
library(tidyverse)

# Number of datasets
n_datasets <- 100

# Set working directory
setwd("/pub59/iasiimwe/TB/datasets")

# Get sex information (Males == 1, Females == 2)
sex <- tibble(fread("/pub59/iasiimwe/TB/Gdata/chr1_filtered.fam", header = FALSE)) %>%
  select(FID = V1, IID = V2) %>%
  left_join(select(tibble(fread("/pub59/iasiimwe/TB/igsr-lwk.tsv.tsv")), IID = `Sample name`, Sex)) %>%
  mutate(Sex = case_when(Sex == "female" ~ 2, Sex == "male" ~ 1))

# Save the sex update file
write.table(sex, file = "sex_update.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Get chromosome 1 (1,255,959 SNPs) as you update sex
plink_command <- "/pub59/iasiimwe/plink2 --bfile /pub59/iasiimwe/TB/Gdata/chr1_filtered --update-sex sex_update.txt --make-bed --out luhya2"
system(plink_command)

# Add folders if they don't exist
if (!dir.exists("true_covar")) dir.create("true_covar")
if (!dir.exists("10_3")) dir.create("10_3")
if (!dir.exists("10_4")) dir.create("10_4")
if (!dir.exists("10_5")) dir.create("10_5")
if (!dir.exists("10_6")) dir.create("10_6")

# Set seed for reproducibility
set.seed(7)

for (i in 1:n_datasets) {
   
  # R2 == 0.1
  # ---------
  # MAFs
  system("/pub59/iasiimwe/plink2 --bfile luhya2 --freq --out luhya")
  
  Sys.sleep(1)
  
  # Randomly select 3 tag SNPs for each of the 3 MAFs (5%, 10%, 20%) 
  tags <- tibble(fread("luhya.afreq")) %>%
    mutate(MAF = round(ALT_FREQS * 100)) %>%
    filter(MAF %in% c(5, 10, 20)) %>%
    select(ID, MAF) 
  
  # Sample one SNP for each MAF (5%, 10%, and 20%)
  sampled_snps <- tags %>% filter(MAF == 5) %>% sample_n(1) %>%
    bind_rows(tags %>% filter(MAF == 10) %>% sample_n(1)) %>%
    bind_rows(tags %>% filter(MAF == 20) %>% sample_n(1))
  
  # Save only the ID column without header and row names
  write.table(sampled_snps$ID, file = "maf_0_1_snps.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # Find the SNPs tagged by these SNPs (at r2 of 0.1) - big ld window to ensure the entire chromosome is covered
  system("/pub59/iasiimwe/plink1.9/plink --bfile luhya2 --r2 --ld-snp-list maf_0_1_snps.txt --ld-window-kb 300000000 --ld-window 999999999 --ld-window-r2 0.1 --out tagged_snps_0_1")
  ## Get tagging SNPs and then exclude them (unique works on adjacent lines so sorting is crucial)
  Sys.sleep(1)
  system("tail -n +2 tagged_snps_0_1.ld | awk '{$1=$1; print}' | cut -f6 -d' ' | grep -F -v -f maf_0_1_snps.txt | sort | uniq > tagged_snps_0_1_unique.txt")
  Sys.sleep(1)
  system("/pub59/iasiimwe/plink1.9/plink --bfile luhya2 --exclude tagged_snps_0_1_exclude.txt --make-bed --out luhya_0_1")
  Sys.sleep(1)
  
  # R2 == 0.5
  # ---------
  # MAFs
  system("/pub59/iasiimwe/plink2 --bfile luhya_0_1 --freq --out luhya")
  Sys.sleep(1)
  
  # Randomly select 3 tag SNPs for each of the 3 MAFs (5%, 10%, 20%) 
  tags <- tibble(fread("luhya.afreq")) %>%
    mutate(MAF = round(ALT_FREQS * 100)) %>%
    filter(MAF %in% c(5, 10, 20)) %>%
    select(ID, MAF) 
  
  # Sample one SNP for each MAF (5%, 10%, and 20%)
  sampled_snps <- tags %>% filter(MAF == 5) %>% sample_n(1) %>%
    bind_rows(tags %>% filter(MAF == 10) %>% sample_n(1)) %>%
    bind_rows(tags %>% filter(MAF == 20) %>% sample_n(1))
  
  # Save only the ID column without header and row names
  write.table(sampled_snps$ID, file = "maf_0_5_snps.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # Find the SNPs tagged by these SNPs (at r2 of 0.5) - big ld window to ensure the entire chromosome is covered
  system("/pub59/iasiimwe/plink1.9/plink --bfile luhya_0_1 --r2 --ld-snp-list maf_0_5_snps.txt --ld-window-kb 300000000 --ld-window 999999999 --ld-window-r2 0.5 --out tagged_snps_0_5")
  Sys.sleep(1)
  ## Get tagging SNPs and then exclude them
  system("tail -n +2 tagged_snps_0_5.ld | awk '{$1=$1; print}' | cut -f6 -d' ' | grep -F -v -f maf_0_5_snps.txt | sort | uniq > tagged_snps_0_5_exclude.txt")
  Sys.sleep(1)
  system("/pub59/iasiimwe/plink1.9/plink --bfile luhya_0_1 --exclude tagged_snps_0_5_exclude.txt --make-bed --out luhya_0_5")
  Sys.sleep(1)
  
  # R2 > 0.5
  # --------
  # MAFs
  system("/pub59/iasiimwe/plink2 --bfile luhya_0_5 --freq --out luhya")
  Sys.sleep(1)
  
  ## Randomly select 3 SNPs for each of the 3 MAFs (5%, 10%, 20%)
  maf_0_5 <- select(tibble(fread("maf_0_5_snps.txt", header = FALSE, sep = " ")), SNP = V1)
  tags <- tibble(fread("luhya.afreq")) %>%
    mutate(MAF = round(ALT_FREQS * 100)) %>%
    filter(MAF %in% c(5, 10, 20)) %>%
    filter(!ID %in% maf_0_5$SNP) %>%
    select(ID, MAF) 
  
  # Sample one SNP for each MAF (5%, 10%, and 20%)
  sampled_snps <- tags %>% filter(MAF == 5) %>% sample_n(1) %>%
    bind_rows(tags %>% filter(MAF == 10) %>% sample_n(1)) %>%
    bind_rows(tags %>% filter(MAF == 20) %>% sample_n(1))
  
  # Save only the ID column without header and row names
  write.table(sampled_snps$ID, file = "maf_0_5_above_snps.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # Find the SNPs tagged by these SNPs (at r2 of 0.5) - big ld window to ensure the entire chromosome is covered
  system("/pub59/iasiimwe/plink1.9/plink --bfile luhya_0_5 --r2 --ld-snp-list maf_0_5_above_snps.txt --ld-window-kb 300000000 --ld-window 999999999 --ld-window-r2 0.5 --out tagged_snps_0_5_above")
  Sys.sleep(1)
  # Get tagging SNPs 
  system("tail -n +2 tagged_snps_0_5_above.ld | awk '{$1=$1; print}' | cut -f6 -d' ' | grep -F -v -f maf_0_5_above_snps.txt | sort | uniq > tagged_snps_0_5_above_include.txt")
  Sys.sleep(1)
  
  # Get remaining SNPs in R
  maf_0_1 <- select(tibble(fread("maf_0_1_snps.txt", header = FALSE, sep = " ")), SNP = V1)
  maf_0_5 <- select(tibble(fread("maf_0_5_snps.txt", header = FALSE, sep = " ")), SNP = V1)
  maf_0_5_above <- select(tibble(fread("maf_0_5_above_snps.txt", header = FALSE, sep = " ")), SNP = V1)
  maf_0_5_above_tagged <- select(tibble(fread("tagged_snps_0_5_above_include.txt", header = FALSE, sep = " ")), SNP = V1)
  maf_0_5_above_ld <- tibble(fread("tagged_snps_0_5_above.ld")) %>%
    filter(!SNP_B %in% maf_0_5_above$SNP) %>%
    select(SNP_A, SNP_B, R2) %>%
    arrange(SNP_A, desc(R2))
  min_rows <- maf_0_5_above_ld %>%
    group_by(SNP_A) %>%
    count() %>% # 155, 5, 9
    pull(n) %>%
    min()
  maf_0_5_above_ld <- maf_0_5_above_ld %>%
    group_by(SNP_A) %>%
    slice_head(n = min_rows) %>%
    pull(SNP_B)
  
  all_include <- c(maf_0_1$SNP, maf_0_5$SNP, maf_0_5_above$SNP, maf_0_5_above_ld)
  length(all_include) # 24
  all_exclude <- c(maf_0_1$SNP, maf_0_5$SNP, maf_0_5_above$SNP, maf_0_5_above_tagged$SNP)
  length(all_exclude) # 178
  true_covar <- c(maf_0_1$SNP, maf_0_5$SNP, maf_0_5_above$SNP)
  length(true_covar) # 9
  
  remaining <- tibble(fread("luhya_0_5.bim")) %>%
    select(SNP = V2) %>%
    filter(!SNP %in% all_exclude)
  nrow(remaining) # 1235843
  
  # To sample 10^3−𝑥, 10^4−𝑥, 10^5−𝑥, and 10^6−𝑥, where 𝑥 = 24, the number of tentatively included SNPs.
  # Sampling function
  sample_rows <- function(dat, num_rows, already_included = 24) { # 24 already included
    # Sample the rows
    sampled_rows <- tibble(SNP = all_include) %>%
      bind_rows(dat[sample(nrow(dat), num_rows - already_included), ]) %>% 
      distinct() 
    
    n_rows <- min(num_rows, nrow(sampled_rows)) # Any cases when more samples are sampled
    sampled_rows <- sampled_rows[1:n_rows, ] 
    
    n_diff <- num_rows - nrow(sampled_rows)
    if (n_diff > 0) { # Few cases when fewer samples are sampled
      unsampled <- dat %>%
        filter(!SNP %in% sampled_rows$SNP)
      sampled_rows <- sampled_rows %>%
        bind_rows(unsampled[1:n_diff, ]) %>% 
        distinct()
    }
    
    while (nrow(sampled_rows) < num_rows) {
      unsampled <- dat %>%
        filter(!SNP %in% sampled_rows$SNP)
      n_to_sample <- num_rows - nrow(sampled_rows)
      sampled_rowsx <- unsampled[sample(nrow(unsampled), n_to_sample), ] 
      sampled_rowsx <- sampled_rowsx[1:n_to_sample, ] # In case more than the required is sampled
      sampled_rows <- sampled_rowsx %>%
        bind_rows(sampled_rows) %>% 
        distinct()
    }
    return(sampled_rows)
  }
  
  # 10^6 dataset
  dat_10_6 <- sample_rows(remaining, num_rows = 10^6) 
  # 10^5 dataset
  dat_10_5 <- sample_rows(dat_10_6, num_rows = 10^5)  
  # 10^4 dataset
  dat_10_4 <- sample_rows(dat_10_5, num_rows = 10^4)
  # 10^3 dataset
  dat_10_3 <- sample_rows(dat_10_4, num_rows = 10^3) 
  
  # Save only the ID column without header and row names
  write.table(dat_10_3$SNP, file = "dat_10_3_snps.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(dat_10_4$SNP, file = "dat_10_4_snps.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(dat_10_5$SNP, file = "dat_10_5_snps.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(dat_10_6$SNP, file = "dat_10_6_snps.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(true_covar, file = paste0("true_covar/true_snps_", i, ".txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # Get the datasets
  system(paste0("/pub59/iasiimwe/plink2 --bfile luhya2 --make-bed --extract dat_10_6_snps.txt --out 10_6/dat_10_6_", i))
  Sys.sleep(1)
  system(paste0("/pub59/iasiimwe/plink2 --bfile 10_6/dat_10_6_", i, " --make-bed --extract dat_10_5_snps.txt --out 10_5/dat_10_5_", i))
  Sys.sleep(1)
  system(paste0("/pub59/iasiimwe/plink2 --bfile 10_5/dat_10_5_", i, " --make-bed --extract dat_10_4_snps.txt --out 10_4/dat_10_4_", i))
  Sys.sleep(1)
  system(paste0("/pub59/iasiimwe/plink2 --bfile 10_4/dat_10_4_", i, " --make-bed --extract dat_10_3_snps.txt --out 10_3/dat_10_3_", i))

  message(paste0("##################\n##################\n##################\n", round(i * 100/n_datasets, 2), "% complete\n##################\n##################\n##################"))
}
# system("wc -l 10_3/dat_10_3_2.bim")



