from sisana.analyze_networks.analyze import calc_tt, calculate_additional_comparison_stats
import sys
import pandas as pd
import pytest

def test_comp():
    # Data for the DataFrame
    data = [[0.4, 0.2, 0.4, 0.1, 0.2, 0.6, 0.7, 0.7, 0.7, 0.7], 
            [0.5, 0.2, 0.4, 0.2, 0.5, 0.6, 0.5, 0.5, 0.5, 0.7], 
            [5, 6, 6, 4, 5, 4, 3, 4, 6, 4], 
            [0.8, 0.8, 1, 1, 0.9, 0.2, 0.3, 0.2, 0.1, 0.3], 
            [1, 0.9, 0.8, 0.9, 1, 0.3, 0.1, 0.3, 0.1, 0.2]]

    # Define row (index) and column names
    genes = ['gene1', 'gene2', 'gene3', 'gene4', 'gene5']
    samps = ['s1t1', 's2t1', 's3t1', 's4t1', 's5t1', 's1t2', 's2t2', 's3t2', 's4t2', 's5t2']
    compdf = pd.DataFrame(data, index=genes, columns=samps)
    
    
    tt_pval = calc_tt(compdf, samps[0:5], samps[5:11], "tt", "tt_pval", "test_statistic")
    tt_res = calculate_additional_comparison_stats(tt_pval, compdf, "group1", "group2", samps[0:5], samps[5:11], "tt", "mediandiff", "tt_pval", "test_statistic")[0]
          
    paired_tt_pval = calc_tt(compdf, samps[0:5], samps[5:11], "paired_tt", "paired_tt_pval", "paired_tt_stat")
    paired_tt_res = calculate_additional_comparison_stats(paired_tt_pval, compdf, "group1", "group2", samps[0:5], samps[5:11], "tt", "mediandiff", "paired_tt_pval", "paired_tt_stat")[0]
    
    wilcox_pval = calc_tt(compdf, samps[0:5], samps[5:11], "wilcoxon", "wilcoxon_pval", "test_statistic")
    wilcox_res = calculate_additional_comparison_stats(wilcox_pval, compdf, "group1", "group2", samps[0:5], samps[5:11], "wilcoxon", "mediandiff", "wilcoxon_pval", "test_statistic")[0]  

    print(wilcox_res.at["gene4", "wilcoxon_pval"])

    assert round(tt_res.at["gene5", "test_statistic"], 3) == 12.348
    assert round(tt_res.at["gene3", "FDR"], 3) == 0.143
    assert round(paired_tt_res.at["gene5", "paired_tt_pval"], 5) == 0.00025
    assert round(paired_tt_res.at["gene2", "median_group2"], 1) == 0.5
    assert round(wilcox_res.at["gene4", "wilcoxon_pval"], 3) == pytest.approx(0.063, abs=0.002) #Note: The real value is something like 0.0624999999999 or something, so we have to approximate here
    assert round(wilcox_res.at["gene2", "median_group1"], 1) == 0.4