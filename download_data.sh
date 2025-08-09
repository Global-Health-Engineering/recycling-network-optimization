#!/bin/bash


echo "Downloading dataset from Zenodo..."
wget -O bin-alloc-data.zip "https://zenodo.org/records/16780698/files/bin-alloc-data.zip?download=1"

# Check if download was successful
if [ $? -eq 0 ]; then
    echo "Download completed successfully."
    
    # Extract to data directory
    echo "Extracting to data/ directory..."
    rm data/raw_data/README.md
    rm data/derived_data/README.md
    unzip -q bin-alloc-data.zip -d data/
    
    # Check if extraction was successful
    if [ $? -eq 0 ]; then
        echo "Extraction completed successfully."
        
        # Clean up the files
        rm bin-alloc-data.zip
        rm -rf data/__MACOSX
        echo "Cleaned up zip file."

    else
        echo "Error: Failed to extract the zip file."
        exit 1
    fi
else
    echo "Error: Failed to download the dataset."
    exit 1
fi

echo "Dataset successfully downloaded and extracted to data/ directory."