#!/bin/bash

# Retrieve cancertype names
cancertypes=($(find "Data Transfer" -mindepth 1 -maxdepth 1 -type d -exec basename {} \;))

# Loop through each cancertype
for cancertype in "${cancertypes[@]}"
do
  # Create RNAseq_cancertype folder
  mkdir -p "RNAseq_${cancertype#*-}"
  
  # Create metadata folder
  mkdir -p "RNAseq_${cancertype#*-}/metadata"
  
  # Extract contents of clinical.cart*.tar.gz to metadata folder
  clinical_tar=$(find "Data Transfer/$cancertype" -name "clinical.cart.*.tar.gz" -type f)
  tar -xzf "$clinical_tar" -C "RNAseq_${cancertype#*-}/metadata" --strip-components=1
  
  # Move gdc_sample_sheet*.tsv to metadata folder
  sample_sheet_file=$(find "Data Transfer/$cancertype" -name "gdc_sample_sheet.*.tsv" -type f)
  if [ -f "$sample_sheet_file" ]; then
    mv "$sample_sheet_file" "RNAseq_${cancertype#*-}/metadata"
  fi
  
  # Create data folder
  mkdir -p "RNAseq_${cancertype#*-}/data"
  
  # Extract files from gdc_download*.tar.gz to data folder
  download_tar=$(find "Data Transfer/$cancertype" -name "gdc_download_*_*.tar.gz" -type f)
  tar -xzf "$download_tar" -C "RNAseq_${cancertype#*-}/data" --strip-components=1
done
