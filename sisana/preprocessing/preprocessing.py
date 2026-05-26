# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import sys
import os
from pathlib import Path
from .file_validation import check_file_extension

def preprocess_data(exp: str, datatype: str, number: int, outdir: str):
    """
    Description:
        This code performs a survival analysis between two user-defined groups and outputs
        both the survival plot and the statistics for the comparison(s)
        
    Parameters:
    -----------
        - exp: str, Path to the expression input data, with row names as genes and column names as samples
        - datatype: str, Type of data (either "csv", "txt", or "tsv")
        - number: int, The number of samples that must express a gene, eitherwise they are removed from downstream analysis
        - outdir: str, Path to output directory
        
    Returns:
    -----------
        - list: [output file path (str), number of genes kept (int), and number of removed genes (int)]
    """
    
    # Create output file prefix by removing the .txt suffix
    expoutfile = Path(exp).stem
    
    # Create output directory if one does not already exist
    os.makedirs(outdir, exist_ok=True)
    
    # Create output for temp files if one does not already exist
    os.makedirs('./tmp/', exist_ok=True)
    
    check_file_extension(exp, datatype)
        
    if datatype == "csv":
        expdf = pd.read_csv(exp, index_col = 0)
    elif datatype == "txt" or datatype == "tsv":
        expdf = pd.read_csv(exp, index_col = 0, sep = "\t")

    with open("./tmp/samples.txt", "w") as f:
        for col in expdf.columns:
            f.write(col + "\n")
      
    # num_originalgenes = len(expdf)
    nsamps = len(expdf.columns)

    num_nonzeros = np.count_nonzero(expdf.iloc[:,0:nsamps], axis=1) # Count the number of non zeros for each gene in the df
    expdf["num_samps_expressed"] = num_nonzeros # add the count to a new col in the df
    
    # Subset df for only genes that appear in at least k samples (user-specified)
    samples_removed = len(expdf[expdf["num_samps_expressed"] < number])
    
    cutdf = expdf[expdf["num_samps_expressed"] >= number]
    cutdf = cutdf.drop(columns=['num_samps_expressed'])

    cutdf.columns = expdf.columns[:-1]
    nsamps_remaining = len(cutdf.columns)    
        
    with open("./tmp/num_samples.txt", "w") as f:
        f.write(str(nsamps_remaining))
    
    # Print summary statistics for the filtering of exp files
    dist = expdf["num_samps_expressed"].value_counts()
    dist = dist.rename_axis('Number of samples with expression').reset_index(name='Number of instances')
    dist = dist.sort_values(by='Number of samples with expression', ascending=False)
    dist['Removed?'] = dist['Number of samples with expression'].apply(lambda x: 'yes' if x < number else 'no')
    print(dist.to_string(index=False))
    sum_kept = dist[dist['Removed?'] == 'no']['Number of instances'].sum()
    sum_removed = dist[dist['Removed?'] == 'yes']['Number of instances'].sum()
    print(f"Number of genes removed in total: {sum_removed}")
    print(f"Number of genes remaining for downstream analysis: {sum_kept}")

    basename = f"{expoutfile}_preprocessed.txt"
    file_outloc = os.path.join(outdir, basename)

    cutdf.to_csv(file_outloc, sep = "\t")    
    print(f"\nFile saved: {file_outloc}")
    
    return(file_outloc, sum_kept, sum_removed)