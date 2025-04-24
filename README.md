# Improving Genotype Imputation and SNP Selection Using Multiple Imputation and Machine Learning 

## Overview
This repository explores and compares **machine learning** and **traditional genome-wide association study (GWAS) approaches** for imputation, covariate selection, and computational efficiency. 

## Background
TB remains a significant global health challenge, causing over 1 million deaths annually. Factors such as genetics, drug resistance, and pharmacokinetics profoundly impact treatment outcomes. The **UNITE4TB consortium** addresses these challenges using advanced tools like **genome-wide association studies (GWAS)** and ML techniques.

Pyrazinamide is integral to TB treatment regimens but shows considerable pharmacokinetic variability due to pharmacogenetic and other influences. This repository uses **simulated pharmacokinetic datasets** to assess the performance of ML and traditional approaches in:
1. Covariate imputation (handling missing SNP data under MCAR, MAR, and MNAR mechanisms).
2. Covariate selection.
3. Computational performance.

## Objectives
The repository aims to:
1. **Compare ML and PMX methods** for high-dimensional pharmacogenetic data:
   - **Covariate Imputation**: Evaluate ML methods (e.g., Random Forest) alongside traditional techniques.
   - **Covariate Selection**: Compare ML approaches (e.g., Penalized regression, ChatGPT) with traditional techniques like GWAS.
   - **Computational Efficiency**: Assess models based on run times, accuracy, bias, precision, and F1 scores.

2. **Provide open-source tools and resources**, including:
   - Scripts for dataset simulation.
   - Methods for covariate imputation and selection.
   - Comparative analysis workflows.

## Dataset Description
The datasets used in this repository are simulated and include:
- **Pharmacokinetic models**: Simulations are based on Vinnard et al.'s one-compartment model with first-order elimination and transit absorption. For reference, see [this publication](https://pmc.ncbi.nlm.nih.gov/articles/PMC5667771/).
- **Genetic data**: SNPs are sourced from the 1000 Genomes Project, incorporating varying effect sizes, allele frequencies, and correlation thresholds.

Key features:
- SNP data with up to 50% missing values introduced under MCAR, MAR, and MNAR mechanisms.
- Nine true covariates among 1,000–1,000,000 SNPs to mimic high-dimensional pharmacogenetic data.

## Analysis Workflow
The analysis is conducted in a **Linux-based command-line environment**, utilizing **R (base)** or **Python (langchain_ai)** environments. Detailed instructions for installing `langchain` are available [here](https://github.com/langchain-ai/langchain). Alternatively, any Python environment can be used. 

### Key Steps
1. Activate your working environment and set the working directories.
2. Ensure all scripts are executable by running: `chmod 755 [script_name]`.
3. Follow these steps:
   ```bash
   # Download Luhya genomic data (only chromosome 1 will be used)
   ./0_download_and_convert.sh

   # Generate binary PLINK datasets containing 1,000, 10,000, 100,000, and 1,000,000 SNPs
   Rscript 1_get_datasets.R

   # Create binary PLINK datasets with missing data (missingness levels: 5%, 10%, 20%, 50%) under MCAR, MAR, MNAR
   Rscript 2_get_missing_datasets.R

   # Generate simulated clinical covariate datasets
   Rscript 3_get_clinical_covariates_dataset.R

   # Simulate concentration-time data
   Rscript 4_get_dv_datasets.R

   # Re-fit base/structural models (without covariates) to obtain EBEs and investigate ETA shrinkage
   Rscript 5_refitting_base_models.R

   # Re-fit models with SNP effects to obtain ‘true’ parameter values (without missing data)
   Rscript 6_refitting_complete_dv_data.R

   # Impute missing data using OpenAI's GPT-4o via API ('ChatGPT' via command line), an imputation server and MICE 
   python 7a_chatgpt.py
   Rscript 7b_Minimac4_preparation_si.R
   Rscript 7b_Minimac4_preparation_mi.R
   Rscript 7c_Minimac4_imputation_server_si.R
   Rscript 7c_Minimac4_imputation_server_mi.R
   ./7d_Minimac4_results_extraction_si.sh
   ./7d_Minimac4_results_extraction_mi_from_si.sh
   ./7d_Minimac4_results_extraction_mi.sh
   ./7e_Minimac4_info_score_extraction_si.sh
   ./7e_Minimac4_info_score_extraction_mi.sh
   Rscript 7f_mice.R

   # Refit imputed datasets using several methods
   Rscript 8a_fit_singly_imputed_data.R 
   Rscript 8b_fit_mice_imputed_data.R
   Rscript 8c_fit_michigan_mi.R
   Rscript 8c_fit_michigan_mi_from_si.R

   # Compare the imputation methods
   Rscript 9_imputation_methods_comparison.R

   # Run covariate selection methods
   Rscript 10_cossac_samba_test.R  
   Rscript 11_gwas_approach.R  
   Rscript 12_ml_methods_1000_with_correlations.R  
   Rscript 12_ml_methods.R

   # Apply dimensionality reduction during covariate selection
   Rscript 13_plink_prune.sh  
   Rscript 13_pruning_optimization.R  
   Rscript 13b_gwas_approach_pruned.R  
   Rscript 13c_ml_methods_pca.R  
   Rscript 13c_ml_methods_pruned_pca.R  
   Rscript 13c_ml_methods_pruned.R

   # Extract key comparison/performance metrics from the results of the covariate selection
   Rscript 14_from_method_tb_1000.R  
   Rscript 14_from_method_tb.R

   # Generate specific figures
   Rscript 15_fig1.R  
   Rscript 15_figs.R  
   Rscript 15_figS15.R

   # Test imputation and covariate selection methods on clinical/diobank datasets
   Rscript 16_iwpc_warfarin.R  
   Rscript 16_ukb.R  
   Rscript 16_warfarin.R
