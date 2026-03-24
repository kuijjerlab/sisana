import pytest
import os
from os.path import join
import subprocess
import yaml
from sisana.preprocessing import preprocess_data
import pandas as pd
from pathlib import Path
from sisana.generate import make_panda_network, make_lioness_networks
from sisana.analyze_networks import compare_bw_groups, survival_analysis, perform_gsea, plot_volcano, plot_expression_degree, plot_heatmap, plot_clustermap, summarize

@pytest.fixture(scope="session")
def shared_temp_dir(tmp_path_factory):
    # Create a base temporary directory for the entire session
    base_dir = tmp_path_factory.mktemp("shared_data")
    # You can add setup code here, like creating shared files
    # (base_dir / "important_file.txt").write_text("This is shared data.")
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
    assert params["generate"]["exp"] == "./output/preprocess/BRCA_TCGA_20_LumA_LumB_samps_5000_genes_exp_preprocessed.txt"
    assert params["generate"]["ppi"] == "./example_inputs/ppi_prior_2024.tsv"
    assert params["generate"]["motif"] == "./example_inputs/motif_prior_names_2024.tsv"
    
def test_preprocess(shared_temp_dir):
    params_path = join(shared_temp_dir, "example_inputs/params.yml")
    params = yaml.load(open(join(params_path)), Loader=yaml.FullLoader)

    os.chdir(shared_temp_dir)
    results = preprocess_data(params["preprocess"]["exp_file"], 
                    params["preprocess"]['filetype'], 
                    params["preprocess"]['number'],
                    params["preprocess"]['outdir'])  
    
    fname, genes_kept, genes_removed = results[0], results[1], results[2] 
    
    assert genes_kept == 4908
    assert genes_removed == 93
   
def test_reconstruct_nw_calc_degree(shared_temp_dir):
    params_path = join(shared_temp_dir, "example_inputs/params.yml")
    params = yaml.load(open(join(params_path)), Loader=yaml.FullLoader)
    os.chdir(shared_temp_dir)
    
    # Make panda networks
    pan = make_panda_network(params["generate"]["exp"],
                                params["generate"]["motif"],
                                params["generate"]["ppi"],
                                params["generate"]["compute"],
                                params["generate"]["modeProcess"],
                                pandafilepath=params["generate"]["pandafilepath"])

    panda_deg_path = f"{str(params['generate']['pandafilepath'])[:-4]}_indegree.csv" # Note: This is the path made in lioness_df_indeg_outdeg_calculator.py
    panda_ind = pd.read_csv(panda_deg_path, nrows=5, index_col=0)
    assert round(panda_ind.loc['A4GNT', 'force'], 2) == -135.86

    # Make lioness networks
    lion_files = make_lioness_networks(panda=pan,
                        compute=params["generate"]["compute"],
                        ncores=params["generate"]["ncores"],
                        start=params["generate"]["start"],
                        end=params["generate"]["end"],
                        lioness_fpath=params["generate"]["lionessfilepath"],
                        panda_fpath=params["generate"]["pandafilepath"])
    
    lioness_ind = pd.read_csv(str(lion_files["lioness_indeg_filepath"]), nrows=5, index_col=0)
    assert round(lioness_ind.loc['A4GNT', 'TCGA_E2_A10E_01A_21R_A10J_07'], 2) == -176.64
    assert round(lioness_ind.loc['AACS', 'TCGA_E2_A14Q_01A_11R_A12D_07'], 2) == -169.59

# def test_compare(shared_temp_dir):
#     params_path = join(shared_temp_dir, "example_inputs/params.yml")
#     params = yaml.load(open(join(params_path)), Loader=yaml.FullLoader)
#     outfiles = compare_bw_groups(datafile=params["compare"]["datafile"], 
#                                 mapfile=params["compare"]["mapfile"], 
#                                 datatype=params["compare"]["datatype"], 
#                                 groups=params["compare"]["groups"],
#                                 testtype=params["compare"]["testtype"], 
#                                 filetype=params["compare"]["filetype"],
#                                 rankby_col=params["compare"]["rankby"],
#                                 outdir=params["compare"]["outdir"])

#     compare_output = pd.read_csv(outfiles[0], nrows=5, index_col=0, sep='\t')
#     assert round(compare_output.loc['RPS4Y2', 'mw_signed_-log(pvalue)'], 2) == 17.33
#     assert round(compare_output.loc['GUCA2B', 'mean_LumB'], 2) == -137.90
