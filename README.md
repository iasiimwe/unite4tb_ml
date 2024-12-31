# Machine Learning vs. Traditional Pharmacometrics in High-Dimensional Pharmacogenetic Pyrazinamide Datasets: Covariate Imputation and Selection

## Overview
This repository explores and compares **machine learning (ML)** and **traditional pharmacometrics (PMX)** approaches for analyzing high-dimensional pharmacogenetic datasets, focusing on **pyrazinamide**, a key drug in tuberculosis (TB) treatment. The primary goals include evaluating covariate imputation, selection methods, and computational efficiency. 

Note: This is an ongoing project. Additional scripts will be uploaded, and existing ones may be optimized as necessary.

## Background
TB remains a significant global health challenge, causing over 1 million deaths annually. Factors such as genetics, drug resistance, and pharmacokinetics profoundly impact treatment outcomes. The **UNITE4TB consortium** addresses these challenges using advanced tools like **genome-wide association studies (GWAS)** and ML techniques.

Pyrazinamide is integral to TB treatment regimens but exhibits considerable pharmacokinetic variability due to pharmacogenetic and other influences. This repository leverages **simulated pharmacokinetic datasets** to assess the performance of ML and PMX approaches in:
1. Covariate imputation (handling missing SNP data under MCAR, MAR, and MNAR mechanisms).
2. Covariate selection.
3. Computational performance.

## Objectives
The repository aims to:
1. **Compare ML and PMX methods** for high-dimensional pharmacogenetic data:
   - **Covariate Imputation**: Evaluate ML methods (e.g., Random Forest, neural networks) alongside traditional techniques.
   - **Covariate Selection**: Compare ML approaches (e.g., penalized regression, autoencoders, ChatGPT) with PMX techniques like Stepwise Covariate Modeling (SCM).
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
- True covariates among 1,000–1,000,000 SNPs to mimic high-dimensional pharmacogenetic data.

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

   # Impute missing data using ChatGPT
   python 7a_chatgpt.py

   # Compare imputation methods (starting with simple ones)
   Rscript 7_imputation_methods_comparison.R
