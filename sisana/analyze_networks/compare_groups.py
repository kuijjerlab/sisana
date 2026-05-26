import os
import pandas as pd
import numpy as np
from pathlib import Path
import csv
from scipy import stats
from scipy.stats import ttest_rel, wilcoxon
from .analyze import file_to_list, map_samples, calc_tt, calculate_additional_comparison_stats
from sisana.exceptions import NotASubsetError
import sys
from numpy import log
from sisana.preprocessing import check_no_hyphens_in_group_names, check_for_header, validate_metadata, check_file_extension

__author__ = 'Nolan Newman'
__contact__ = 'nolankn@uio.no'
    
def compare_bw_groups(datafile: str, mapfile: str, datatype: str, groups: list, testtype: str, data_filetype: str, map_filetype: str, rankby_col: str, outdir: str):
    '''
    Description:
        This code compares the degree or expression of two sample groups
     
    Parameters:
    -----------
        - datafile: str, Path to the data file
        - mapfile: str, Path to the mapping file, which maps sample name to sample group
        - datatype: str, The type of data being used ("expression" or "degree")
        - groups: str, Names of the two groups (from the second column of mapfile) to be compared. The second group listed will be used as the numerator in the fold change calulation.
        - testtype: str, Type of comparison to perform, either "tt" for Student's t-test, "mw" for Mann-Whitney U, "paired_tt", or "wilcoxon"
        - data_filetype: str, The type of data file ("csv" or "txt" or "tsv") being used
        - map_filetype: str, The type of mapping file ("csv" or "txt" or "tsv") being used
        - rankby_col: str, Choices: ["mediandiff", "mwu", "neglogp", "meandiff"]. The statistic to rank the .rnk output file by for GSEA. 
        - outdir: str, The directory to save the output to
        
    Returns:
    -----------
        - list of the output file paths
    '''
    
    # Create output directory if one does not already exist    
    os.makedirs(outdir, exist_ok=True)
    
    check_file_extension(datafile, data_filetype)
    check_file_extension(mapfile, map_filetype)
    
    if data_filetype == "csv":
        datadf = pd.read_csv(datafile, index_col = 0)
    else:
        datadf = pd.read_csv(datafile, index_col = 0, sep = "\t")
        
    if map_filetype == "csv":
        mapfile = pd.read_csv(mapfile, index_col = 0)
    else:
        mapfile = pd.read_csv(mapfile, index_col = 0, sep = "\t")
    
    check_for_header(datafile, data_filetype)
    validate_metadata(mapfile)
        
    if testtype == "tt" or testtype == "mw":
        check_no_hyphens_in_group_names(mapfile)
        
        # Assign samples from mapping file to groups
        sampdict = map_samples(mapfile, groups[0], groups[1])
        total_samps = len(sampdict[groups[0]]) + len(sampdict[groups[1]])
    
    elif testtype == "paired_tt" or testtype == "wilcoxon":
        sampdict = {}
        sampdict[groups[0]] = list(mapfile.index)
        sampdict[groups[1]] = mapfile.iloc[:,0].tolist()
        total_samps = len(sampdict[groups[0]]) + len(sampdict[groups[1]])
                   
    # remove unnecessary samples to save on memory
    if len(datadf.columns) > total_samps:
        allsamps = []
        allsamps = sampdict[groups[0]] + sampdict[groups[1]]
        compdf = datadf.loc[:, allsamps]
    else:
        compdf = datadf
        
    del datadf
    
    # Validate that the samples in the mapping file are a subset of those in the data frame
    samps_list_group1 = sampdict[groups[0]]
    samps_list_group2 = sampdict[groups[1]]
    
    print(f"Numbers of samples in group 1: {len(samps_list_group1)}")
    print(f"Numbers of samples in group 2: {len(samps_list_group2)}")
    
    if not set(samps_list_group1).issubset(list(compdf.columns)):
        raise NotASubsetError(user_list=samps_list_group1, data_list=compdf.columns, dtype="samples")
    if not set(samps_list_group2).issubset(list(compdf.columns)):
        raise NotASubsetError(user_list=samps_list_group2, data_list=compdf.columns, dtype="samples")

    print("Performing comparisons, please wait...")
    
    # Calculate p-value/FDR
    
    pval_column = testtype + "_pvalue"
    test_stat_column = "test_statistic"

    if testtype != "mw": 
        newpvaldf = calc_tt(compdf, sampdict[groups[1]], sampdict[groups[0]], testtype, pval_column, test_stat_column)
        
    else:
        # mwu_calculations = calc_tt(compdf, sampdict[groups[1]], sampdict[groups[0]], testtype, pval_column, test_stat_column)
        mwu_calculations = compdf.apply(lambda row : calc_tt(compdf, row[sampdict[groups[1]]], row[sampdict[groups[0]]], testtype, pval_column, test_stat_column), axis = 1)
    
        # Format the output data frame
        pval_column = testtype + "_pvalue"
        
        test_stat_column = "mw_uvalue"
        neglogp_column = "mw_signed_-log(pvalue)"
        cles_column = "mw_CLES"
        
        pvaldf = pd.DataFrame({'Target':mwu_calculations.index, pval_column:mwu_calculations})
        newpvaldf = pd.DataFrame(pvaldf[pval_column].to_list(), columns=[test_stat_column, pval_column, neglogp_column, cles_column])
        newpvaldf['Target'] = mwu_calculations.index
        newpvaldf = newpvaldf.set_index('Target')

    # For expression, do regular log2FC
    # For degrees, first need to transform the edge value by doing ln(e^w + 1),
    # then calculate degrees. Then you can do the log2FC of degrees
    # Calculate log2FC (only for expression data, since you can have negative 
    # degree values). This transformation is described in the paper "Regulatory Network 
    # of PD1 Signaling Is Associated with Prognosis in Glioblastoma Multiforme"
    # if datatype == "expression":
    #     fc = compdf.apply(lambda row : calc_log2_fc(row[sampdict[groups[1]]], row[sampdict[groups[0]]]), axis = 1)

    # Calculate additional statistics (difference of means, difference of medians, etc.) and add to the output data frame
    if testtype != "mw":
        newpvaldf, ranked = calculate_additional_comparison_stats(newpvaldf, compdf, groups[0], groups[1], samps_list_group1, samps_list_group2, testtype, rankby_col, pval_column, test_stat_column)
    else:
        newpvaldf, ranked = calculate_additional_comparison_stats(newpvaldf, compdf, groups[0], groups[1], samps_list_group1, samps_list_group2, testtype, rankby_col, pval_column, test_stat_column, neglogp_column=neglogp_column, cles_column=cles_column)
    # Write to disk
    if testtype == "tt" or testtype == "mw":
        save_file_path = os.path.join(outdir, f"comparison_{testtype}_between_{groups[0]}_{groups[1]}_{datatype}.txt")
        
        if testtype == "tt":
            save_file_path_ranked = os.path.join(outdir, f"comparison_{testtype}_between_{groups[0]}_{groups[1]}_{datatype}_ranked_test_stat.rnk")
        else:
            save_file_path_ranked = os.path.join(outdir, f"comparison_{testtype}_between_{groups[0]}_{groups[1]}_{datatype}_ranked_{rankby_col}.rnk")
            
    if testtype == "paired_tt" or testtype == "wilcoxon":
        save_file_path = os.path.join(outdir, f"comparison_{testtype}_{datatype}.txt")
        save_file_path_ranked = os.path.join(outdir, f"comparison_{testtype}_{datatype}_ranked_{rankby_col}.rnk")
    
    # Make output directory if it does not already exist
    Path(outdir).mkdir(parents=True, exist_ok=True)    
    
    # Remove the abs() columns (previously for sorting) for exporting the file
    columns_to_keep = [x for x in newpvaldf.columns if x[0:3] != "abs"]
    
    newpvaldf.to_csv(save_file_path, columns = columns_to_keep, sep = "\t")
    ranked.to_csv(save_file_path_ranked, sep = "\t", header = False)

    print(f"\nFile saved: {save_file_path}\nThis file contains all the statistics results and is just for your reference.\n")
    print(f"File saved: {save_file_path_ranked}\nThis file contains only the calculated test statistic for each gene and is used for input to the GSEA analysis step.\n")
    
    return([save_file_path, save_file_path_ranked])
