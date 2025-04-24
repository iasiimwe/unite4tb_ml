# Load required packages
# ---------------------
library(tidyverse)
library(data.table)
library(mice) 

select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate

# Add folders if they don't exist
if (!dir.exists("mice")) dir.create("mice")

# Additive coding function - to impute as a continuous covariate, which is consistent with IMPUTE
additive_coding_fn <- function(x) {
  x <- gsub(" ", "", x) 
  x[x == "00"] <- NA
  y <- table(x)
  y <- names(sort(-y)) # Sort in descending order, assume mutant type is the least abundant
  if (length(y) == 3) x <- gsub(y[[3]], "2", x) # Give mutant homozygotes "2"
  if (length(y) > 1) x <- gsub(y[[2]], "1", x) # Give heterozygotes "1" 
  x <- gsub(y[[1]], "0", x) # Give homozygotes 0
  return(as.numeric(x))
}

# MAF formula
maf_fn <- function(x) {
  x[x > 2] <- 2 # In case we have anything rounded to above 0
  x <- table(x)
  if (length(x) == 1) maf <- 1 
  if (length(x) == 2) maf <- as.numeric(x["1"] / (2 * ( x["0"] +  x["1"])))
  if (length(x) == 3) maf <- as.numeric((x["1"] + 2 *x["2"])/ (2 * ( x["0"] +  x["1"] +  x["2"])))
  return(round(maf * 100))
}

# Wrapper function (to help skip to the next iteration if mice fails)

mice_with_skip <- function(data, m, method, pred, seed = 7, printFlag = FALSE, maxit = 20, post = NULL) {
  result <- tryCatch({
    mice(data, m = m, method = method, pred = pred, seed = seed, 
         printFlag = printFlag, maxit = maxit, post = post)
  }, error = function(e) {
    cat("Error in mice_with_skip:", conditionMessage(e), "\n")
    NULL  # Return NULL on error
  })
  # Return the result (NULL if an error occurred)
  return(result)
}

n_datasets <- 100

# Method	    Description	                            Scale Type
# pmm	        Predictive mean matching	              Any 
# midastouch	Weighted predictive mean matching	      Any
# sample	    Random sample from observed values	    Any
# cart	      Classification and regression trees	    Any
# rf	        Random forest imputation	              Any
# mean	      Unconditional mean imputation	          Numeric
# norm	      Bayesian linear regression	            Numeric
# ri	        Random indicator for nonignorable data	Numeric

imputation_methods <- c("pmm", "midastouch", "sample", "cart", "rf", "mean", "norm", "ri")
mechanisms <- c("MCAR", "MAR", "MNAR")
missing_percentages <- c(5, 10, 20, 50)
to_post <- FALSE

for (k in seq_along(imputation_methods)) {
  imputation_method <- imputation_methods[[k]]
  
  for (mechanism in mechanisms) {
    for (missing_percentage in missing_percentages) {
      message(paste0("Starting\n     Method: ",
                     imputation_method, "\n     Mechanism: ", 
                     mechanism, "\n     missing %: ", missing_percentage))
      
      # Times storage matrices
      time_tb <- tibble(end = vector("character", n_datasets), 
                        start = vector("character", n_datasets))
      if (mechanism == "MNAR" && to_post) {
       time_tb_post <- tibble(end = vector("character", n_datasets), 
                              start = vector("character", n_datasets))
      }
      
      for (k in 1:n_datasets) {
        # Get missing SNP data
        dat_path <- paste0("missing/", missing_percentage, "/dat_", mechanism, "_", k)
        
        # Convert to ped format
        system(paste0("/pub59/iasiimwe/plink1.9/plink --bfile ", dat_path, " --recode tab --out temp_", mechanism))
        
        # Process ped and map_files
        ped <- tibble(fread(paste0("temp_", mechanism, ".ped"), header = FALSE))
        map <- tibble(fread(paste0("temp_", mechanism, ".map"), header = FALSE))
        snps <- pull(map, V2)
        ped_start <- ped %>% select(V1:V6) # This selects FID, IID, father's ID, mother's ID, sex and phenotype data
        snp_columns <- colnames(ped)[!colnames(ped) %in% colnames(ped_start)]
        ped_snps <- ped %>%
          select(all_of(snp_columns)) 
        colnames(ped_snps) <- snps 
        
        # True covariates
        true_snps <- fread(paste0("true_covar/true_snps_", k, ".txt"), header = FALSE, sep = " ")$V1
        
        # Bring the true SNPs to the first columns and change to 'additive' coding
        ped_snps <- ped_snps %>%
          select(any_of(true_snps), everything()) %>%
          mutate_all(additive_coding_fn) %>%
          mutate(sex = ped_start$V5, .before = 1)
        
        # Change the SNP names as these can be problematic (we already know the true are the first 9)
        colnames(ped_snps) <- c("sex", paste0("SNP", 1:1000))
        
        # Get time of code execution start
        start <- Sys.time()
        
        # Impute the NAs here
        #---------------------
        # SNPs to always include - in practice, we need all predictors, so we could predict in batches, by looping inlist
        # Batches could be informed by LD - assuming collinearity will not be an issue
        # Note that GWAS's cover one SNP at a go, so as long as the other covariates are included here, the analysis model will be a subset of the imputation model
        inlist <- c("sex", paste0("SNP", 1:9)) # covariates that should be included in every imputation model

        # Let us randomly select 20 other SNPs before applying quick pred
          # The most important thing you can do is downsize the number of predictors in each model to the most important, say, 20-30 variables (https://stefvanbuuren.name/fimd/sec-modelform.html#sec:predictors).
        set.seed(7)
        to_include <- paste0("SNP", sample(c(10:1000), 20))
        ped_snps <- ped_snps %>%
          select(all_of(c(inlist, to_include)))
        
        # Quick selection of predictors from the data (in practice, we need all predictors, so we could predict in batches, by looping inlist)
        pred <- quickpred(ped_snps, include = inlist)
        
        # Path to save
        path_to_save <- paste0("mice/", imputation_method, "_", missing_percentage, "_", 
                               mechanism, "_", k)
        
        imp <- mice_with_skip(ped_snps, m = missing_percentage, method = imputation_method,
                    pred = pred, maxit = 0)
        if (is.null(imp)) next
        
        # Apply post-processing rules in a loop for SNP1 to SNP9 (round to the nearest whole number)
        for (snp in paste0("SNP", 1:9)) {
          imp$post[[snp]] <- "imp[[j]] <- round(imp[[j]])"
        }
        imp <- mice_with_skip(ped_snps, m = missing_percentage, method = imputation_method,
                    pred = pred, maxit = 20, post = imp$post)
        if (is.null(imp)) next
        # stripplot(imp, sex + SNP1 + SNP2 + SNP3 + SNP4 + SNP5 + SNP6 + SNP7 + SNP8 + SNP9 ~ .imp)
        # plot(imp, paste0("SNP", 1:3))
        # complete(imp, 1)
        write_rds(imp, paste0(path_to_save, ".rds")) # ID orders unchanged so no need to add them now
        
        end <- Sys.time()
        
        if (mechanism == "MNAR" && to_post) {
          start_post <- Sys.time()
          inlist <- c("sex", paste0("SNP", 1:9)) # covariates that should be included in every imputation model
          set.seed(7)
          to_include <- paste0("SNP", sample(c(10:1000), 20))
          ped_snps <- ped_snps %>%
            select(all_of(c(inlist, to_include)))
          pred <- quickpred(ped_snps, include = inlist)
          imp <- mice_with_skip(ped_snps, m = missing_percentage, method = imputation_method,
                      pred = pred, maxit = 0)
          if (is.null(imp)) next
          
          # Let us find delta, i.e. a value to add to the imputed values
          delta <- c(0, 0.5, 0.75, 0.8) # 0 is MAR
          imp.all.undamped <- vector("list", length(delta))
          
          sample_d <- function(x) sample(c(0, 1), size = 100, replace = TRUE, prob = c(1-x, x))
          
          for (i in seq_along(delta)) {
            d <- delta[i]
            cmd <- paste0("imp[[j]][,i] <- round(imp[[j]][,i] + sample_d(", d, "))[1:length(imp[[j]][,i])]")
            for (snp in paste0("SNP", 1:9)) {
              imp$post[[snp]] <- cmd
            }
            imp <- mice_with_skip(ped_snps, m = 5, method = imputation_method, # Five imputations for the sensitivity analysis, only want MAF
                        pred = pred, maxit = 20, post = imp$post)
            if (is.null(imp)) next
            imp.all.undamped[[i]] <- imp
            message(paste0("Delta search: ", i, " out of ", length(delta),  " complete!"))
          }
          
          delta_vct <- vector("double", length(delta))
          for (i in seq_along(delta)) {
            if (is.null(imp.all.undamped[[i]])) {
              delta_vct[i] <- 1000000 # A very high value for failed iterations
            } else {
              delta_vct[i] <- imp.all.undamped[[i]] %>%  
                complete("long") %>%
                select(SNP1:SNP9) %>%
                mutate_all(maf_fn) %>%
                distinct() %>%
                mutate_at(vars(SNP1, SNP4, SNP7), function(x) abs(x - 5)) %>%
                mutate_at(vars(SNP2, SNP5, SNP8), function(x) abs(x - 10)) %>%
                mutate_at(vars(SNP3, SNP6, SNP9), function(x) abs(x - 20)) %>%
                rowMeans()
            }
          }
          
          # Use the delta that gives us the values closest to the known MAFs (we know these MAFs or could have used an imputed dataset)
          d <- delta[which.min(delta_vct)]
          cmd <- paste0("imp[[j]][,i] <- round(imp[[j]][,i] + sample_d(", d, "))[1:length(imp[[j]][,i])]")
          for (snp in paste0("SNP", 1:9)) {
            imp$post[[snp]] <- cmd
          }
          
          imp <- mice_with_skip(ped_snps, m = missing_percentage, method = imputation_method,
                      pred = pred, seed = 7, printFlag = FALSE, maxit = 20, post = imp$post)
          if (is.null(imp)) next
          write_rds(imp, paste0(path_to_save, "_post.rds"))
          end_post <- Sys.time()
        } 
        
        time_tb$start[k] <- paste("T", as.character(ymd_hms(start)))
        time_tb$end[k] <- paste("T", as.character(ymd_hms(end)))
        
        if (mechanism == "MNAR" && to_post) {
          time_tb_post$start[k] <- paste("T", as.character(ymd_hms(start_post)))
          time_tb_post$end[k] <- paste("T", as.character(ymd_hms(end_post)))
        }
        
        message(paste0("##################\n##################\n##################\n", 
                       round(k * 100/n_datasets, 2), 
                       "% complete\n     Method: ",imputation_method, "\n     Mechanism: ", 
                       mechanism, "\n     missing %: ", missing_percentage, 
                       "\n##################\n##################\n##################"))
      }
      path_to_save <- paste0("mice/", imputation_method, "_", missing_percentage, "_", mechanism)
      write.csv(time_tb, paste0(path_to_save, "_time.csv"), row.names = FALSE)
      if (mechanism == "MNAR" && to_post) write.csv(time_tb_post, paste0(path_to_save, "_time_post.csv"), row.names = FALSE)
    }
  }
}
