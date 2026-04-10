import pytest
import os
from os.path import join
import subprocess
import yaml
from sisana.preprocessing import preprocess_data
import pandas as pd
from pathlib import Path
import sisana
import sys
from sisana.reconstruct import make_panda_network, make_lioness_networks
from sisana.analyze_networks import compare_bw_groups, survival_analysis, perform_gsea, plot_volcano, plot_expression_degree, plot_heatmap, plot_clustermap, summarize

@pytest.fixture(scope="session")
# Create a base temporary directory for the entire session
def shared_temp_dir(tmp_path_factory):
    base_dir = tmp_path_factory.mktemp("shared_data")
    return base_dir

@pytest.fixture
def get_zenodo_files(shared_temp_dir):
    subprocess.run(["sisana", "-e"], cwd=shared_temp_dir)
    print(f"tmp path: {shared_temp_dir}")
    return join(shared_temp_dir, "example_inputs")

def test_fetch(get_zenodo_files):
    assert os.path.exists(join(get_zenodo_files, 'params.yml'))

# @pytest.fixture
# def load_params(tmp_path):
#     params_path = join(get_zenodo_files, "params.yml")
#     params = yaml.load(open(join(tmp_path, "example_inputs", "params.yml")), Loader=yaml.FullLoader)
#     print(params)
#     return params

def test_load_params(shared_temp_dir):
    params_path = join(shared_temp_dir, "example_inputs/params.yml")
    params = yaml.load(open(join(params_path)), Loader=yaml.FullLoader)
    
    assert params["preprocess"]["exp_file"] == "./example_inputs/BRCA_TCGA_20_LumA_LumB_samps_5000_genes_exp.tsv"    
    assert params["reconstruct"]["exp"] == "./output/preprocess/BRCA_TCGA_20_LumA_LumB_samps_5000_genes_exp_preprocessed.txt"
    assert params["reconstruct"]["ppi"] == "./example_inputs/ppi_prior_2024.tsv"
    assert params["reconstruct"]["motif"] == "./example_inputs/motif_prior_names_2024.tsv"
    
def test_preprocess(shared_temp_dir):
    subprocess.run(["sisana",  "preprocess", "./example_inputs/params.yml"], cwd=shared_temp_dir)
    log_files_path = shared_temp_dir / "log_files"
    print(os.path.join(log_files_path, "preprocess_log.txt"))
    log_yaml = yaml.load(open(os.path.join(log_files_path, "preprocess_log.txt")), Loader=yaml.FullLoader)
    
    # fname, genes_kept, genes_removed = results[0], results[1], results[2] 
    assert log_yaml["Additional information"]["genes kept"] == 4908
    assert log_yaml["Additional information"]["genes removed"] == 93
   
def test_reconstruct_nw_calc_degree(shared_temp_dir):
    subprocess.run(["sisana",  "reconstruct", "./example_inputs/params.yml"], cwd=shared_temp_dir)
    output_files_path = shared_temp_dir / "output/network/"
    panda_ind = pd.read_csv(os.path.join(output_files_path, "panda_network_indegree.csv"), nrows=5, index_col=0)
    assert round(panda_ind.loc['A4GNT', 'force'], 2) == -135.86

    lioness_ind = pd.read_csv(os.path.join(output_files_path, "lioness_indegree.csv"), nrows=5, index_col=0)
    assert round(lioness_ind.loc['A4GNT', 'TCGA_E2_A10E_01A_21R_A10J_07'], 2) == -176.64
    assert round(lioness_ind.loc['AACS', 'TCGA_E2_A14Q_01A_11R_A12D_07'], 2) == -169.59

def test_compare(shared_temp_dir):
    subprocess.run(["sisana",  "compare", "./example_inputs/params.yml"], cwd=shared_temp_dir)
    output_files_path = shared_temp_dir / "output/compare_means/"
    compare_output = pd.read_csv(os.path.join(output_files_path, "comparison_mw_between_LumA_LumB_degree.txt"), nrows=5, index_col=0, sep="\t")
   
    assert round(compare_output.loc['RPS4Y2', 'mw_signed_-log(pvalue)'], 2) == 17.33
    assert round(compare_output.loc['GUCA2B', 'mean_LumB'], 2) == -137.90
