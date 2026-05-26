import scipy 
from scipy import stats
from scipy.stats import ttest_rel, wilcoxon
import scipy.stats
import csv
import re
import pandas as pd
from statistics import mean, median
import math
import numpy as np
import sys
from sisana.exceptions import NotASubsetError

def file_to_list(fname):
    """
    This function takes as input a text file with one object per line and returns the contents of that file as a list

    Args:
        fname: the name of the file to convert to a list

    Returns:
        returnlist: list of objects from the text file
    """
    all_lines = open(fname, "r").read().splitlines()
    returnlist = [name for name in all_lines if name]

    return (returnlist)   

def map_samples(mapfile: pd.DataFrame, type1: str, type2: str):
    '''
    Function that assigns samples to groups for statistical analysis

        Arguments:
            - mapfile: pd.DataFrame, data frame with sample name as rows, first column contains the name of the group each sample belongs to
            - type1: str, the name of the group in the first set of samples, must be present in the mapfile
            - type2: str, the name of the group in the second set of samples, must be present in the mapfile
        
        Returns:
            - samp_type_dict: dict, A dictionary that contains group name as keys and a list of sample names that belong to that group as values 
    '''

    # Check if the supplied groups are a subset of the mapping column
    # mapf = pd.read_csv(mapfile, index_col = 0)
    mapf = mapfile
    input_groups_set = set([type1, type2])
    column_group_set = set(mapf[mapf.columns[0]])

    if not input_groups_set.issubset(column_group_set):
        raise NotASubsetError([type1, type2], mapf[mapf.columns[0]], "groups")

    samp_type_dict = {group: mapf.index[mapf.iloc[:,0] == group].tolist() for group in mapf.iloc[:,0].unique()}
    return (samp_type_dict)

def calc_tt(compdf, group1, group2, ttype, pcol, tcol):
    '''
    Performs either a students t-test or mann-whitney test between two groups
    
        Arguments:
            - compdf: data frame containing values for comparisons
            - group1: list of samples in first group
            - group2: list of samples in second group
            - ttype: str, the type of test to perform, either tt or mw
            - pcol: str, the name of the column containing the p-value in the output df
            - tcol: str, the name of the column containing the test statistic in the output df
                 
    '''

    if ttype == 'tt':
        pval = compdf.apply(lambda row : stats.ttest_ind(row[group1], row[group2]), axis = 1)
        
        pvaldf = pd.DataFrame({'Target':pval.index, pcol:pval.values})
        newpvaldf = pd.DataFrame(pvaldf[pcol].to_list(), columns=['test_statistic', pcol])
        newpvaldf['Target'] = pval.index
        newpvaldf = newpvaldf.set_index('Target')
        return(newpvaldf)

    elif ttype == 'paired_tt':
        tval = []
        pval = []

        for gene, row in compdf.iterrows():
            x = row[group1].values.astype(float)
            y = row[group2].values.astype(float)

            t, p = ttest_rel(x, y)
            tval.append(t)
            pval.append(p)

        newpvaldf = pd.DataFrame({'Target':compdf.index, pcol:pval, tcol: tval})
        newpvaldf.set_index('Target', inplace=True)
        return(newpvaldf)

    elif ttype == 'wilcoxon':
        tval = []
        pval = []
        
        for gene, row in compdf.iterrows():
            x = row[group1].astype(float).values
            y = row[group2].astype(float).values
            
            t, p = wilcoxon(x, y, method="auto", correction = False)
            tval.append(t)
            pval.append(p)
            
        newpvaldf = pd.DataFrame({'Target':compdf.index, pcol:pval, tcol: tval})
        newpvaldf.set_index('Target', inplace=True)
        return(newpvaldf)
        
    elif ttype == 'mw':        
        from pingouin import mwu

        mw_results = mwu(group1, group2)
        mw_results = mw_results.iloc[0, :].values.tolist()
        mw_results_flipped = mwu(group2, group1)
        mw_results_flipped = mw_results_flipped.iloc[0, :].values.tolist()
        statistic_group1, pval, cles = mw_results[0], mw_results[2], mw_results[4]
        statistic_group2 = mw_results_flipped[0]
        
        mwu_min_ustat = min(statistic_group1, statistic_group2)

        # Calculate the -log(pval) for ranking of genes and GSEA
        neglog_pval = -1 * math.log(pval)

        if median(group1) < median(group2):
            neglog_pval = neglog_pval * -1            

        return(mwu_min_ustat, pval, neglog_pval, cles)  

    

    
def calculate_additional_comparison_stats(newpvaldf, compdf, name_group1, name_group2, samps_group1, samps_group2, testtype, rankby_col, pval_column, test_stat_column, *args, **kwargs):
    '''
    Calculates the difference of means across two groups
        
        Arguments:
            - group1: list of samples in first group
            - group2: list of samples in second group
            - difftype: 
    '''
    def _calc_group_difference(group1, group2, difftype=["mean", "median"]):
        '''
        Calculates the difference of means across two groups
            
            Arguments:
                - group1: list of samples in first group
                - group2: list of samples in second group
                - difftype: 
        '''
        if difftype == "mean":
            return(mean(group2) - mean(group1))
        elif difftype == "median":
            return(median(group2) - median(group1))
        
    mean_diff = compdf.apply(lambda row : _calc_group_difference(row[samps_group1], row[samps_group2], difftype="mean"), axis = 1)
    median_diff = compdf.apply(lambda row : _calc_group_difference(row[samps_group1], row[samps_group2], difftype="median"), axis = 1)

    print("Comparisons finished...") 
    
    # Calcuate means per group
    mean_g2_colname = f"mean_{name_group2}"    
    mean_g1_colname = f"mean_{name_group1}"      
    newpvaldf[mean_g2_colname] = compdf[samps_group2].mean(axis=1)
    newpvaldf[mean_g1_colname] = compdf[samps_group1].mean(axis=1)
    meandiff_colname = f"difference_of_means_({name_group2}-{name_group1})"      
    newpvaldf[meandiff_colname] = mean_diff
    newpvaldf["abs(difference_of_means)"] = abs(mean_diff)

    # Calcuate medians per group    
    median_g2_colname = f"median_{name_group2}"    
    median_g1_colname = f"median_{name_group1}"      
    newpvaldf[median_g2_colname] = compdf[samps_group2].median(axis=1)
    newpvaldf[median_g1_colname] = compdf[samps_group1].median(axis=1)
    mediandiff_colname = f"difference_of_medians_({name_group2}-{name_group1})"      
    newpvaldf[mediandiff_colname] = median_diff
    newpvaldf["abs(difference_of_medians)"] = abs(median_diff)

    # Perform multiple test correction
    FDR_colname = "FDR"
    newpvaldf[FDR_colname] = stats.false_discovery_control(newpvaldf[pval_column])
    newpvaldf = newpvaldf.sort_values(pval_column, ascending = True)
    
    if testtype == "mw": 
        if rankby_col == "mwu":
            sortcol = "mw_uvalue"
        elif rankby_col == "mediandiff":
            sortcol = f"difference_of_medians_({name_group2}-{name_group1})"
        elif rankby_col == "meandiff":
            sortcol = f"difference_of_means_({name_group2}-{name_group1})"
        elif rankby_col == "neglogp":
            sortcol = "mw_signed_-log(pvalue)"
    else:
        sortcol = test_stat_column
    
    # Create new df without pval, ranked on test statistic (as chosen by user)
    ranked = newpvaldf.sort_values(sortcol, ascending = False)
    ranked.drop([pval_column, FDR_colname], inplace=True, axis=1)
    ranked = ranked[sortcol]
    
    # Rearrange column order so that FDR calculations comes after p-value
    if testtype != "mw":
        colorder = [test_stat_column, pval_column, FDR_colname,
                    mean_g2_colname, mean_g1_colname,
                    meandiff_colname, median_g2_colname, median_g1_colname,
                    mediandiff_colname]
        newpvaldf = newpvaldf.loc[:, colorder] 
    else:
        neglogp_column = kwargs.get('neglogp_column')
        cles_column = kwargs.get('cles_column')
        colorder = [test_stat_column, pval_column, neglogp_column, FDR_colname,
                    cles_column, mean_g2_colname, mean_g1_colname,
                    meandiff_colname, median_g2_colname, median_g1_colname,
                    mediandiff_colname]
        newpvaldf = newpvaldf.loc[:, colorder]
        
    return([newpvaldf, ranked])

############################################################################
### Note: The following functions were removed after determing that 
###       calculating the log fold change in edges is not a good metric,
###       since the edges are already logged to begin with. May rethink
###       this calculation and re-add them later
############################################################################
# def transform_edge_to_positive_val(edgeval):
#     '''
#     Transform degree values to be positive so fold change can be calculated. This transformation 
#     is described in the paper "Regulatory Network of PD1 Signaling Is Associated with Prognosis 
#     in Glioblastoma Multiforme"
        
#         Arguments:
#             - edgeval: float, value of edge
            
#         Returns:
#             - float, transformed value of edge
#     '''
#     # We get an inf if we have too large of values, but since transforming a large value
#     # does not change the resulting value (i.e. ln(e^1000) + 1) = 1000, we can avoid the
#     # inf values by just not transforming those edges
#     if edgeval > 700:
#         newval = edgeval
#     else:
#         newval = math.log(np.exp(edgeval) + 1)
#     if newval == float("inf"):
#         raise Exception(f"val {edgeval} gives a result of inf")
#     return(newval)

# def calc_log2_fc(group1, group2):
#     '''
#     Calculates the log2 fold change of means across two groups
        
#         Arguments:
#             - group1: list of samples in first group
#             - group2: list of samples in second group
#     '''
    
#     # print("\n")
#     # print(f"group1 values:")
#     # print(group1)
#     # print(f"Mean of group 1: {scipy.mean(group1)}")
#     # print(f"group2 values:")
#     # print(group2)
#     # print(f"Mean of group 2: {scipy.mean(group2)}")

#     ### Note: For the following calculations, assuming the user has followed the previous SiSaNA steps,
#     ### neither group will ever have a mean of 0 since the preprocess.py step filters out any genes of 
#     ### low abundance. So to simplify the logic for this step, I have removed any checks for group means
#     ### of 0. May need to implement this change later though if people are not running preprocess.py.

#     # if scipy.mean(group1) == 0 and scipy.mean(group2) == 0:
#     #     log2FC = 0
#     # elif scipy.mean(group1) == 0:
        
#     # elif scipy.mean(group2) != 0:
    
#     avg_g1 = mean(group1)
#     avg_g2 = mean(group2)

#     if avg_g1 == 0 or avg_g2 == 0:
#         log2FC = "NA"
#         # print("log2FC is NA")

#     elif all(i >= 0 for i in group1) and all(i >= 0 for i in group2):
#         log2FC = math.log2(avg_g2/avg_g1)
#         # print(f"log2FC is {log2FC}\n")

#         # print(log2FC)
#     else:
#         print(group1, group2)
#         raise Exception("\n\nError: Negative values found in data. The log2 fold change can only be calculated on expression data. Negative values indicate that degree was likely used as an input type instead.\n")
#     # elif scipy.mean(group2) == 0:
#     #     # Take the smallest non-zero value, divide by 10, and use that as the average value for the denominator
#     #     no_zeros = [i for i in group1 if i != 0]
#     #     denom = min(no_zeros)/10
#     #     log2FC = math.log2(scipy.mean(group1)/denom)        
#     # else:
#     #     log2FC = "NA"
    
#     return(log2FC)
