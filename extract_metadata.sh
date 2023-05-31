#!/bin/bash

# Retrieve cancertype names
cancertypes=($(find "Data Transfer" -mindepth 1 -maxdepth 1 -type d -exec basename {} \;))

# Loop through each cancertype
for cancertype in "${cancertypes[@]}"
do
  # Extract clinical.cart*.tar.gz to temporary folder
  extracted_folder=$(find "Data Transfer/$cancertype" -name "clinical.cart.*.tar.gz" -type f)
  echo "Extracted folder: $extracted_folder"  # Debugging info
  tmp_folder="RNAseq_${cancertype#*-}/metadata/tmp"
  mkdir -p "$tmp_folder"
  tar -xzf "$extracted_folder" -C "$tmp_folder" --strip-components=1
  
  # Move extracted files to metadata folder
  echo "Files to move:"
  find "$tmp_folder" -type f -print  # Debugging info
  find "$tmp_folder" -type f -exec mv -t "RNAseq_${cancertype#*-}/metadata/" {} +
  
  # Clean up temporary folder
  rm -r "$tmp_folder"
  
  # Remove the original compressed folder
  rm "$extracted_folder"
done
