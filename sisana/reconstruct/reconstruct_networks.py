
import numpy as np
import pandas as pd
import argparse
from netZooPy.panda.panda import Panda
from netZooPy.lioness.lioness import Lioness
from pathlib import Path
from sisana.postprocessing import convert_lion_to_pickle
from sisana.analyze_networks import calculate_panda_degree, calculate_lioness_degree

import os

def make_panda_network(exp: str, motif: str, ppi: str, compute: str, modeProcess: str, pandafilepath: str):    
    '''
    Reconstructs the PANDA network, using netZooPy      
     
    Parameters:
    -----------
        - exp: str, Path to preprocessed expression file
        - motif: str, Path to motif file
        - ppi: str, Path to PPI file
        - compute: str, Whether to use cpu or gpu
        - modeProcess: str, From netZooPy documentation: "legacy" refers to the processing mode in netZooPy<=0.5. "union" n of all TFs and genes across priors and fills the missing genes in the priors with zeros. "intersection"  input genes and TFs across priors and removes the missing TFs/genes.
        - pandafilepath: str, Path to the PANDA output file to be created. Must have a .txt extension.
        
    Returns:
    -----------
        - PANDA object
    '''


    pandapath = Path(pandafilepath)
    if str(pandapath)[-4:] != ".txt":
        raise Exception("Error: Panda output file must have a .txt extension. Please edit your pandafilepath variable in your params file.")
    os.makedirs(pandapath.parent, exist_ok=True)
        
    panda_obj = Panda(exp, 
            motif_file=motif, 
            ppi_file=ppi, 
            computing=compute,
            modeProcess=modeProcess,
            save_tmp=False, 
            remove_missing=False, 
            keep_expression_matrix=True, 
            save_memory=False,
            with_header=True)
    
    panda_res = panda_obj.export_panda_results
    panda_res.to_csv(pandafilepath, sep=" ", index=False)
            
    print("Now calculating PANDA degrees...")
    calculate_panda_degree(inputfile=pandafilepath)
    
    return(panda_obj)

def make_lioness_networks(panda, compute: str, ncores: str, start, end, lioness_fpath: str, panda_fpath: str):    
    '''
    Reconstructs LIONESS networks, using netZooPy        
     
    Parameters:
    -----------
        - panda: object returned by export_panda_results() when running make_panda_network()
        - compute: str, Whether to use "gpu" or "cpu" for computing networks
        - ncores: str, Number of cores to use for calculating LIONESS networks. Must be less than or equal to the number of samples.
        - start: To start from nth sample. The background is always what is used for PANDA and stays the same. If this value and the end value are not specified, then all samples will be used.
        - end: To end at nth sample. Must be used when a value for "start" is specified. The background is always what is used for PANDA and stays the same
        - lioness_fpath: str, Path to the lioness output file to be created. Must have a .npy extension
        - panda_fpath: str, Path to panda file reconstructed with make_panda_network()
        
    Returns:
    -----------
        - PANDA object
    '''
    if lioness_fpath[-4:] != ".npy":
        raise Exception("Error: Lioness output file must have a .npy extension. Please edit your lionessfilepath variable in your params file.")

    lionesspath_no_ext = lioness_fpath[:-4]

    # If user wants to run lioness in batches or only run for certain samples (e.g. 10 samples at a time), then do the following.
    # Note that this still uses all samples as the background, but will only reconstruct networks for the given sample numbers
    if start is not None:
        lionesspath_new_path = Path(f"{lionesspath_no_ext}_samples_{start}_to_{end}.npy")                
    else:
        lionesspath_new_path = Path(lioness_fpath)

    lioness_full_path = Path(lioness_fpath)

    # Run Lioness on a subset of samples if specified in the params file, otherwise run on all samples
    if start is not None:
        Lioness(panda, 
                computing=compute, 
                precision="double",
                ncores=ncores, 
                save_dir=lioness_full_path.parent, 
                save_fmt="npy",
                start=start,
                end=end)
    else:
        Lioness(panda, 
                computing=compute, 
                precision="double",
                ncores=ncores, 
                save_dir=lioness_full_path.parent, 
                save_fmt="npy")
        
    os.rename(os.path.join(lioness_full_path.parent, "lioness.npy"), lionesspath_new_path)
    
    # Convert lioness numpy file to pickle file
    print("\nLIONESS networks created. Now converting results to a .pickle file...")

    liondf = pd.DataFrame(np.load(lionesspath_new_path))   
    
    if start is not None:
        pickle_path = f"./tmp/lioness_samples_{start}_to_{end}.pickle"
    else:
        pickle_path = './tmp/lioness.pickle'
    
    liondf = convert_lion_to_pickle(panda=panda_fpath,
                    lion=liondf,
                    type="npy", 
                    names='./tmp/samples.txt',  
                    outfile=pickle_path,
                    start=start,
                    end=end)
    
    # Calculate degrees
    print("\n.pickle file created. Now calculating LIONESS degrees...")
    calculate_lioness_degree(nwdf=liondf,
                                pickle=pickle_path)
    print("LIONESS degrees have now been calculated.")

    if start is not None:
        lioness_indeg_filename = f"lioness_indegree_samples_{start}_to_{end}"
        lioness_indeg_filename_path = f"{Path(lioness_full_path).parent}/{lioness_indeg_filename}.csv"
        lioness_outdeg_filename = f"lioness_outdegree_samples_{start}_to_{end}"
        lioness_outdeg_filename_path = f"{Path(lioness_full_path).parent}/{lioness_outdeg_filename}.csv"

        Path(f"./tmp/lioness_samples_{start}_to_{end}_indegree.csv").rename(lioness_indeg_filename_path)
        Path(f"./tmp/lioness_samples_{start}_to_{end}_outdegree.csv").rename(lioness_outdeg_filename_path)

    else:
        lioness_indeg_filename = f"lioness_indegree"
        lioness_indeg_filename_path = f"{Path(lioness_full_path).parent}/{lioness_indeg_filename}.csv"
        lioness_outdeg_filename = f"lioness_outdegree"
        lioness_outdeg_filename_path = f"{Path(lioness_full_path).parent}/{lioness_outdeg_filename}.csv"

        Path("./tmp/lioness_indegree.csv").rename(lioness_indeg_filename_path)
        Path("./tmp/lioness_outdegree.csv").rename(lioness_outdeg_filename_path)

    print(f"LIONESS network saved to {str(lionesspath_new_path)}")
    print(f"LIONESS degrees saved to:")
    print(f"{Path(lioness_full_path).parent}/{lioness_indeg_filename}.csv")
    print(f"{Path(lioness_full_path).parent}/{lioness_outdeg_filename}.csv")

    lion_outputs = {}
    lion_outputs["lioness_nw_filepath"]  = lionesspath_new_path
    lion_outputs["lioness_indeg_filepath"]  = lioness_indeg_filename_path
    lion_outputs["lioness_outdeg_filepath"]  = lioness_outdeg_filename_path
    return(lion_outputs)

    