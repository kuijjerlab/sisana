import numpy as np
import pandas as pd
import argparse
from .post import files_to_dfs
import sys
from sisana.postprocessing import convert_lion_to_pickle
from sisana.analyze_networks import calculate_lioness_degree
import qnorm
from pathlib import Path

def quantile_normalize_edges(net: str, pandafilepath: str, outdir: str, start, end):    
    '''
    Quantile normalizes the edges of the lioness network and saves the output in a pickle file
     
    Parameters:
    -----------
        - net: str, lioness npy file created in the "reconstruct" step
        - outdir: str, Path to output directory
        
    Returns:
    -----------
        - Nothing
    '''

    liondf = np.load(net)
    
    print("Quantile normalizing edges, please wait...")
    normalized_data = pd.DataFrame(qnorm.quantile_normalize(liondf, axis=1))
    
    print("Quantile normalization complete, now calculating degree...")

    # Note: need to convert to pickle file first here as its required for use in calculate_lioness_degree(). Not ideal for speed,
    # since we need to read in the pickle file afterwards too, but sufficient for now and can be changed later.

    pickle_outpath = f"{Path(net).parent}/{Path(net).stem}_quantile_normalized.pickle"
    file_path_and_stem = pickle_outpath.split(".")[0]
    outdeg_filename =  f"{file_path_and_stem}_outdegree.csv"
    indeg_filename =  f"{file_path_and_stem}_indegree.csv"

    liondf = convert_lion_to_pickle(panda=pandafilepath,
                lion=normalized_data,
                type="npy", 
                names='./tmp/samples.txt',  
                outfile=pickle_outpath,
                start=start,
                end=end)
    
    calculate_lioness_degree(nwdf=liondf,
                                pickle=pickle_outpath)
    print("LIONESS degrees have now been calculated.")
    print(f"Normalized LIONESS network saved to {str(pickle_outpath)}")
    print(f"Indegrees calculated from the normalized LIONESS network saved to {str(indeg_filename)}")
    print(f"Outdegrees calculated from the normalized LIONESS network saved to {str(outdeg_filename)}")

    sys.exit(0) 