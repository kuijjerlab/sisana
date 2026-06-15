import numpy as np
import pandas as pd
import argparse
from .post import files_to_dfs
import sys
from sisana.postprocessing import convert_lion_to_pickle
from sisana.analyze_networks import calculate_lioness_degree
import qnorm
from pathlib import Path

def quantile_normalize_edges(net: str, filetype: str, pandafilepath: str, num_cpus: int, start, end):    
    '''
    Quantile normalizes the edges of the lioness network and saves the output in a pickle file
     
    Parameters:
    -----------
        - net: str, lioness npy file created in the "reconstruct" step
        - filetype: str, file type of lioness file, either npy or h5
        - num_cpus: int, number of CPUs to use for quantile normalization (4 is generally the max recommended)
        
    Returns:
    -----------
        - Nothing
    '''

    print(f"Loading LIONESS network from {net}")
    if filetype == "npy":
        liondf = np.load(net)
    else: 
        liondf = pd.read_hdf(net, key="Basal")
            
    print("Quantile normalizing edges, please wait...")
    normalized_data = pd.DataFrame(qnorm.quantile_normalize(liondf, axis=1, ncpus=num_cpus))
    
    print("Quantile normalization complete, now calculating degree...")

    # Note: need to convert to pickle file first here as its required for use in calculate_lioness_degree(). Not ideal for speed,
    # since we need to read in the pickle file afterwards too, but sufficient for now and can be changed later.

    pickle_outpath = f"{Path(net).parent}/{Path(net).stem}_quantile_normalized_edges.pickle"
    file_path_and_stem = pickle_outpath.split(".")[0]
    outdeg_filename =  f"{file_path_and_stem}_outdegree.csv"
    indeg_filename =  f"{file_path_and_stem}_indegree.csv"
    
    if filetype == "npy":
        print(f"Saving quantile normalized LIONESS network to {pickle_outpath} for use in degree calculation...")
        liondf = convert_lion_to_pickle(panda=pandafilepath,
                    lion=normalized_data,
                    type="npy", 
                    names='./tmp/samples.txt',  
                    outfile=pickle_outpath,
                    start=start,
                    end=end)
    else: # Assuming that the h5 file already has the sample names as column names and TF-gene edges as row names, so we can skip the convert_lion_to_pickle step and just save the quantile normalized data frame as a pickle file for use in degree calculation.
        normalized_data.to_pickle(pickle_outpath)

    print("Calculating LIONESS degrees from the quantile normalized network...")
    
    if filetype == "npy":
         calculate_lioness_degree(nwdf=liondf,
                                  pickle=pickle_outpath)
    else:
        calculate_lioness_degree(nwdf=normalized_data,
                                 pickle=pickle_outpath)
        
    print("LIONESS degrees have now been calculated.")
    print(f"\nNormalized LIONESS network saved to {str(pickle_outpath)}")
    print(f"Indegrees calculated from the normalized LIONESS network saved to {str(indeg_filename)}")
    print(f"Outdegrees calculated from the normalized LIONESS network saved to {str(outdeg_filename)}")
    
    return([pickle_outpath, indeg_filename, outdeg_filename])
