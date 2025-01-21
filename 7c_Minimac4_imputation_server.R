# # To log in: open command prompt
# ssh iasiimwe@euler01.liv.ac.uk
# ssh gauss06
# source /pub59/iasiimwe/miniconda3/bin/activate base
# cd /pub59/iasiimwe/imputationserver
# export PATH=/pub59/iasiimwe/imputationserver/singularity/go/bin:$PATH
# export PATH=/pub59/iasiimwe/imputationserver/singularity/bin:$PATH
# export PATH=/pub59/iasiimwe/imputationserver/docker:$PATH
# export PATH=/pub59/iasiimwe/imputationserver:$PATH
# export PATH=/pub59/iasiimwe/imputationserver/jdk-23.0.1+11/bin:$PATH
# # The local web service can now be started. By default, it runs on port 8082.
# java -Xmx10G -Djavax.net.ssl.trustStore=/pub59/iasiimwe/imputationserver/jdk-23.0.1+11/lib/security/cacerts -jar cloudgene.jar server
# # In a new terminal on your local machine, create an SSH tunnel through euler01 to gauss06
# ssh -L 8082:gauss06:8082 iasiimwe@euler01.liv.ac.uk
# # Once the tunnel is established, open a browser and go to:
# http://localhost:8082

# Load libraries
library(data.table)
library(tidyverse) 

# Relevant function
not_scientific <- function(x) trimws(format(x, scientific = FALSE))

TOKEN <- "eyJhbGciOiJIUzI1NiJ9.eyJzdWIiOiJhZG1pbiIsIm5iZiI6MTczNzIwMjMxNSwibWFpbCI6ImkuYXNpaW13ZUBsaXZlcnBvb2wuYWMudWsiLCJhcGlfaGFzaCI6InRxVm9FS2tFUVpPc01XVHExcDFsNUlOUnhQRGdpciIsInJvbGVzIjpbXSwiaXNzIjoiY2xvdWRnZW5lIiwibmFtZSI6ImFkbWluIiwiYXBpIjp0cnVlLCJleHAiOjM4ODQ2ODU5NjIsInRva2VuX3R5cGUiOiJBUElfVE9LRU4iLCJpYXQiOjE3MzcyMDIzMTUsInVzZXJuYW1lIjoiYWRtaW4ifQ.VkVCVJtWgkwdaNl11f0d1Qmfv3OvIHBSeoYpiqwPUj4"

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

