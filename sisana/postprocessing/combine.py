import os
import pandas as pd
import numpy as np
from pathlib import Path
import csv
import sys
import glob
from sisana.postprocessing import convert_lion_to_pickle

__author__ = 'Nolan Newman'
__contact__ = 'nolankn@uio.no'
    
def combine_files(params_dict: str):
    '''
    Description:
        This code combines the batched indegree and outdegree files into one file for 
        each, and optionally combines the network files as well. It also deletes the intermediate 
        files if specified in the params.yml file.
     
    Parameters:
    -----------
        - params_dict: str, the path to the updated params.yml file

    Returns:
    -----------
        - list of the output file paths
    '''

    degree_dir_path = str(Path(params_dict['degree_dir']))
    
    with open('./tmp/samples.txt', 'r') as file:
        samplist = file.readlines()
        samplist = [samp.strip() for samp in samplist] 

    panda_file = pd.read_csv(params_dict['panda_file'], sep = " ")
    panda_file["edge"] = panda_file["tf"] + "-" + panda_file["gene"] 
    panda_file.index = panda_file["edge"]

    def _get_batched_files(regex: str, ext: str):
        """
        Description:
            Finds the batched indegree and outdegree files, saving them to their own lists

        Parameters:
        -----------     
            - regex: str, regular expression to use for finding files
            - ext: str, the extension of the files 
        
        Returns:
        -----------
            - List of the created files
        """
        df_list = []
        filenames_list = []
        print("Files found to combine:")
        if ext == "csv":
            for file in glob.glob(f"{degree_dir_path}/{regex}"):
                print(f"  - {file}")
                df = pd.read_csv(file, index_col=0)
                df_list.append(df)
                filenames_list.append(file)
        else:
            for file in glob.glob(f"{degree_dir_path}/{regex}"):
                print(f"  - {file}")
                noext = file[:-4]
                startsamp, endsamp = int(noext.split("_")[-3]), int(noext.split("_")[-1])

                numpy_file = np.load(file)
                data = pd.DataFrame(numpy_file)

                data.columns = samplist[startsamp-1:endsamp]
                data.index = panda_file.index
        
                df_list.append(data)
                filenames_list.append(file)
                
        return(df_list, filenames_list)

    # Combine the degree files automatically, since they are relatively small.
    # Combining networks may run into memory issues, so it's optional
    print(f"\nCombining indegree files, please wait...")
    indeg_dataframes, indeg_filenames = _get_batched_files("lioness_indegree_samples_*_to_*.csv", "csv")
    combined_indeg = pd.concat(indeg_dataframes, axis=1)
    indeg_filepath = f"{degree_dir_path}/lioness_indegree.csv"
    combined_indeg.to_csv(indeg_filepath, index=True)
    print(f"File created: {indeg_filepath}")
    
    print(f"\nCombining outdegree files, please wait...")
    outdeg_dataframes, outdeg_filenames = _get_batched_files("lioness_outdegree_samples_*_to_*.csv", "csv")       
    combined_outdeg = pd.concat(outdeg_dataframes, axis=1)
    outdeg_filepath = f"{degree_dir_path}/lioness_outdegree.csv"
    combined_outdeg.to_csv(outdeg_filepath, index=True)
    print(f"File created: {outdeg_filepath}")
                
    if params_dict["delete_intermediate_files"] == True:
        [os.remove(file) for file in indeg_filenames]
        [os.remove(file) for file in outdeg_filenames]
                
    if params_dict["networks"] == True:   
        print(f"\nCombining network files, please wait...")
        numpy_dataframes, numpy_filenames = _get_batched_files("lioness_networks_samples_*_to_*.npy", "npy")                  
        combined_nw = pd.concat(numpy_dataframes, axis=1)
                    
        pickle_path = './tmp/lioness.pickle'
        np_path = f"{degree_dir_path}/lioness_network.npy"
        
        with open("./tmp/combined_samples.txt", "w") as f:
            for col in combined_nw.columns:
                f.write(col + "\n")

        convert_lion_to_pickle(params_dict['panda_file'],
                    combined_nw,
                    "npy", 
                    './tmp/combined_samples.txt',  
                    pickle_path)
        
        combined_nw = combined_nw.to_numpy()
        np.save(np_path, combined_nw)
        
        # combined_nw.to_csv(np_path, index=True)
        print(f"File created: {np_path}")
            
        if params_dict["delete_intermediate_files"] == True:
            [os.remove(file) for file in numpy_filenames]
            
    return([indeg_filepath, outdeg_filepath, np_path])
                    