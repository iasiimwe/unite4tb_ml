# cd /pub59/iasiimwe/TB/datasets

import re 
import os
import io
import pandas as pd
import numpy as np
import time 
import requests
from datetime import datetime
from pathlib import Path
from openai import OpenAI
from PIL import Image
from io import BytesIO
from config import set_environment
set_environment() # Set environment variables from the configuration file

# Dominant coding function
def dominant_coding_fn(x):
    x = x.str.replace(" ", "", regex=False)
    x = x.replace("00", np.nan)
    y = x.value_counts().sort_values(ascending=False).index.tolist()
    if len(y) == 3:
        x = x.replace(y[2], "1")
    if len(y) > 1:
        x = x.replace(y[1], "1")
    x = x.replace(y[0], "0")
    return pd.to_numeric(x, errors="coerce")

# Prompt
prompt = """
Task Description:
You are an expert data scientist specializing in data analysis and imputation. Your task is to analyze and process a tuberculosis (TB) dataset provided as a CSV file. The dataset contains the following columns:

- **Participant ID**: A unique identifier for 100 participants.
- **Sex**: Encoded as 1 for male and 2 for female.
- **Single Nucleotide Polymorphisms (SNPs)**: Columns labeled SNP1 to SNP1000. Among these, the first nine SNPs (SNP1-SNP9) are the covariates of primary interest.

### Missing Data
Some SNP data contain missing values, denoted as NaN. The missing data follow three distinct mechanisms:
1. **Missing Completely at Random (MCAR)**: Missingness is unrelated to the data.
2. **Missing at Random (MAR)**: Missingness is related to other observed data.
3. **Missing Not at Random (MNAR)**: Missingness depends on the value of the missing data itself.

### Objectives
1. **Identify Missingness Mechanisms**:
   - Determine the missingness mechanism (MCAR, MAR, or MNAR) for each SNP1-SNP9.

2. **Impute Missing Values**:
   - Impute the missing values for SNP1-SNP9 using methods that align with the identified missingness mechanisms.
   - Ensure all imputed values are either:
     - `0` (wild-type homozygotes)
     - `1` (mutant allele).

3. **Export Results**:
   - Display the imputed dataset with the following columns: ID, sex, and SNP1-SNP9.
   - Include all 100 participants in your output.
   - SNP1-SNP9 should not have any missing entries.

Example:
ID      sex   SNP1 SNP2 SNP3 SNP4 SNP5 SNP6 SNP7 SNP8 SNP9
NA19017 2     0    0    0    0    0    1    0    0    0
NA19020 1     0    0    1    0    0    0    0    0    0
...

4. **Report Findings**:
   - Provide a summary of the identified missingness mechanisms for SNP1-SNP9.
   - Clearly outline how each imputation method was selected and implemented.
   - Include a link to download the imputed dataset in the format:
     `[Download the imputed dataset](sandbox:/mnt/data/imputed_tb_data.csv)`.

### Requirements
- Imputation methods must strictly correspond to the identified missingness mechanisms.
- Generate clean and reproducible code with explanations in comments.
- Ensure the output aligns with the given specifications.
"""

# Parameters
n_datasets = 100
imputation_methods = ["ChatGPT"]
mechanisms = ["MCAR", "MAR", "MNAR"]
missing_percentages = [5, 10, 20, 50]

# Main process
for imputation_method in imputation_methods:
    for mechanism in mechanisms:
        for missing_percentage in missing_percentages:
            print(f"Starting\n     Method: {imputation_method}\n     Mechanism: {mechanism}\n     Missing %: {missing_percentage}")
            
            # Path to save
            path_to_save = Path(f"chatgpt/{imputation_method}_{missing_percentage}_{mechanism}")
            # path_to_save.mkdir(parents=True, exist_ok=True)
            
            # Results storage
            time_tb = pd.DataFrame({"start": [None] * n_datasets, "end": [None] * n_datasets})
            
            for j in range(1, n_datasets + 1):
                # Get missing SNP data
                dat_path = f"/pub59/iasiimwe/TB/datasets/missing/{missing_percentage}/dat_{mechanism}_{j}"
                
                # Convert to ped format using PLINK
                os.system(f"/pub59/iasiimwe/plink1.9/plink --bfile {dat_path} --recode tab --out temp")
                
                # Process ped and map files
                ped = pd.read_csv("temp.ped", sep="\t", header=None)
                map_df = pd.read_csv("temp.map", sep="\t", header=None)
                snps = map_df.iloc[:, 1].tolist()
                ped_start = ped.iloc[:, :6]  # FID, IID, father's ID, mother's ID, sex, phenotype
                ped_snps = ped.iloc[:, 6:]
                ped_snps.columns = snps
                
                # True covariates
                true_snps = pd.read_csv(f"/pub59/iasiimwe/TB/datasets/true_covar/true_snps_{j}.txt", header=None, sep=" ").iloc[:, 0].tolist()
                
                # Apply 'dominant' coding function to all columns
                ped_snps = ped_snps.apply(dominant_coding_fn, axis=0)
                
                # Remove special characters from true_snps
                true_snps = [re.sub(r"[^a-zA-Z0-9 ]", "", snp) for snp in true_snps]
                
                # Remove special characters and "true" from column names of ped_snps
                ped_snps.columns = [re.sub(r"[^a-zA-Z0-9 ]|true", "", col) for col in ped_snps.columns]
                
                # Select true SNPs first, then retain other columns, preserving the order
                ped_snps = ped_snps[[col for col in true_snps if col in ped_snps.columns] + [col for col in ped_snps.columns if col not in true_snps]]
                
                # Rename columns to SNP1, SNP2, ..., SNP1000
                ped_snps.columns = [f"SNP{i+1}" for i in range(len(ped_snps.columns))]
                
                # Add the 'sex' column to ped_snps
                ped_snps["sex"] = ped_start.iloc[:, 4] 
                ped_snps["ID"] = ped_start.iloc[:, 0] 
                
                # Relocate "ID" and 'sex' column to the first position
                columns = ["ID", "sex"] + [col for col in ped_snps.columns if col not in ["ID", "sex"]]
                ped_snps = ped_snps[columns]
                
                # Save this for upload by assistant
                ped_snps.to_csv("temp_ped_snps.csv", index=False)
                
                # Imputation
                start_time = datetime.now()
                
                print("Connecting to ChatGPT.")
                
                # Initialize client
                client = OpenAI()
                
                # Upload a file with an "assistants" purpose
                file = client.files.create(file=open("temp_ped_snps.csv", "rb"), purpose='assistants')
                
                # Create an assistant using the file ID
                assistant = client.beta.assistants.create(
                    instructions=prompt,
                    model="gpt-4o",
                    temperature=0,
                    tools=[{"type": "code_interpreter"}],
                    tool_resources={"code_interpreter": {"file_ids": [file.id]}}
                )
                
                # Create a thread and add a message for analysis
                thread = client.beta.threads.create(
                    messages=[
                        {
                            "role": "user",
                            "content": "Start analyzing the TB dataset and provide the results.",
                            "attachments": [{"file_id": file.id, "tools": [{"type": "code_interpreter"}]}]
                        }
                    ]
                )
                
                # Loop until a dataset can be extracted
                print("Beginning loop to identify extracted file id.")
                while True:
                    # Stream thread run
                    with client.beta.threads.runs.stream(
                        thread_id=thread.id,
                        assistant_id=assistant.id,
                    ) as stream:
                        stream.until_done()
                    
                    # Obtain messages and extract file paths
                    messages = client.beta.threads.messages.list(thread_id=thread.id)
                    sorted_messages = sorted(messages.data, key=lambda msg: msg.created_at)
                    
                    file_id_extracted = None  # Initialize file ID as None
                    
                    # Iterate through the sorted messages
                    for message in sorted_messages:
                        # Iterate through the content blocks in each message
                        for content_block in message.content:
                            if hasattr(content_block, 'text') and hasattr(content_block.text, 'annotations'):
                                # Look for file path annotations in the text
                                for annotation in content_block.text.annotations:
                                    if annotation.type == "file_path":  # Check if it's a file path annotation
                                        file_id_extracted = annotation.file_path.file_id  # Extract file ID
                    # If a file ID was extracted
                    if file_id_extracted:
                        print(f"Extracted file ID: {file_id_extracted}")
                        
                        # Download the file using the OpenAI client
                        file_data = client.files.content(file_id_extracted)
                        csv_data = file_data.read().decode("utf-8")
                        df = pd.read_csv(io.StringIO(csv_data))
                        
                        # Check for required columns
                        required_columns = ["ID", "sex"] + [f"SNP{i}" for i in range(1, 10)]
                        if all(col in df.columns for col in required_columns):
                            # Identify SNP columns
                            snp_columns = [col for col in df.columns if col.startswith("SNP")]
                            
                            # Check for missing or infinite values in SNP columns
                            if not df[snp_columns].isna().any().any() and not df[snp_columns].isin([float("inf"), -float("inf")]).any().any():
                                # Convert SNP columns to integers
                                try:
                                    df[snp_columns] = df[snp_columns].astype(int)
                                    
                                    # Save the dataset
                                    df.to_csv(f"{path_to_save}_{j}.csv", index=False)
                                    print(f"Dataset saved successfully at {path_to_save}_{j}.csv")
                                    break  # Exit the loop if everything is successful
                                    
                                except Exception as e:
                                    print(f"Error converting to integers: {e}")
                            else:
                                print("Dataset contains missing or infinite values. Retrying...")
                        else:
                            print("Dataset does not contain required columns. Retrying...")
                    else:
                        print("No file ID found. Retrying...")
                    
                    continue  # Repeat the loop if conditions are not met  
                
                # List to store images
                images = []
                
                # Loop through all sorted messages
                for message in sorted_messages:
                    # Loop through the content blocks in each message
                    for content_block in message.content:
                        if content_block.type == 'image_file':  # Check if it's an image
                            file_id = content_block.image_file.file_id  # Extract file ID
                            
                            # Download the image using the file ID
                            image_data = client.files.content(file_id)
                            image = Image.open(BytesIO(image_data.read()))  # Open image
                            images.append(image)  # Add to the list
                
                # Combine images into a single file
                if images:
                    # Calculate total height and max width for the combined image
                    total_height = sum(image.height for image in images)
                    max_width = max(image.width for image in images)
                    
                    # Create a blank image with the required dimensions
                    combined_image = Image.new("RGB", (max_width, total_height))
                    
                    # Paste each image into the combined image
                    y_offset = 0
                    for img in images:
                        combined_image.paste(img, (0, y_offset))
                        y_offset += img.height
                        
                    # Save the combined image as a .png file
                    combined_image.save(f"{path_to_save}_{j}_combined_analysis_images.png")
                    print("Combined image saved as 'combined_analysis_images.png'")
                else:
                    print("No images found in the messages.")
                
                # Save messages to a log file
                output_file = f"{path_to_save}_{j}_log.txt"    
                with open(output_file, "w") as log_file:
                    for message in sorted_messages:
                        for content_block in message.content:
                            if hasattr(content_block, 'text') and hasattr(content_block.text, 'value'):
                                text_value = content_block.text.value
                                log_file.write(text_value + "\n")  
                
                print(f"Messages saved to {output_file}")
                
                # Log end time and save data
                end_time = datetime.now()                
                time_tb.loc[j - 1, "start"] = start_time.strftime("%Y-%m-%d %H:%M:%S")
                time_tb.loc[j - 1, "end"] = end_time.strftime("%Y-%m-%d %H:%M:%S")
                
                print(f"{round(j * 100 / n_datasets, 2)}% complete\n     Method: {imputation_method}\n     Mechanism: {mechanism}\n     Missing %: {missing_percentage}")
            
            # Save time log
            time_tb.to_csv(f"{path_to_save}_time.csv", index=False)
            time.sleep(3)
