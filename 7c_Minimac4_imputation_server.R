# Load libraries
library(data.table)
library(tidyverse) 

# Relevant function
not_scientific <- function(x) trimws(format(x, scientific = FALSE))

TOKEN <- "abc123..."

n_datasets <- 100

for (ref_population in c("afr", "eur")) {
  # ID storage tibble (time will be tracked from the logs)
  id_tb <- tibble(dataset = c(1:n_datasets), 
                  job_id = vector("character", n_datasets))
  
  for (k in 1:n_datasets) {
    # Path
    vcf_file_path <- paste0("/pub59/iasiimwe/TB/datasets/michigan/chr", k, ".vcf.gz")
    
    # Construct the curl command
    curl_command <- paste0(
      'curl http://localhost:8082/api/v2/jobs/submit/imputationserver2 ',
      '-H "X-Auth-Token: ', TOKEN, '" ',
      '-F "files=@', vcf_file_path, '" ',
      '-F "refpanel=1000g-phase-3-v5-public" ',
      '-F "population=', ref_population, '" ',
      '-F "password=abc123" ',
      '-F "mode=imputation" ',
      '-F "phasing=eagle" ',
      '-F "build=hg19" ',
      '-F "r2Filter=0"'
    )
    
    # Capture the output of the curl command
    response <- system(curl_command, intern = TRUE)
    
    # Parse the response to extract the job id
    response_json <- jsonlite::fromJSON(paste(response, collapse = ""))
    
    # Extract the job id
    id_tb$job_id[k] <- response_json$id 
    
    message(paste0("##################\n##################\n##################\n", 
                   toupper(ref_population), ": ", round(k * 100/n_datasets, 2), "% complete!", 
                   "\n##################\n##################\n##################"))
    
    # Check if k is a multiple of 5 (we want to run at most 5 jobs on the imputation server)
    if (k %% 5 == 0) {
      message("k is a multiple of 5. Pausing for 8 minutes...")
      Sys.sleep(8 * 60) # Pause for 8 minutes (test runs took about 6.5 minutes i.e. wait for 5 runs to complete before submitting again)
      message("Resuming execution.")
    } 
  }
  write.csv(id_tb, paste0("/pub59/iasiimwe/TB/datasets/michigan/job_id_", ref_population, ".csv"), row.names = FALSE)
} 

