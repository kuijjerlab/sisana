
import numpy as np
import pandas as pd
import argparse

def convert_lion_to_pickle(panda: str, lion: str, type: str, names: str, outfile: str, start: int=None, end: int=None, hotstart: bool=False):    
    '''
    Creates data frames from the input panda and lioness files        
     
    Parameters:
    -----------
        - panda: str, Path to panda output file
        - lion: pd.DataFrame, lioness data frame
        - type: str, file type of lioness file, either npy or txt
        - names: str, File with list of sample names (one per line) in the same order that were supplied for panda/lioness
        - outfile: str, Path to output file in pickle format (e.g. lioness.pickle)
        - start: the starting index of samples that networks were created for. Only specify this value if "start" and "stop" params were used in the params.yml file.
        - end: the ending index of samples that networks were created for. Only specify this value if "start" and "stop" params were used in the params.yml file.
        
    Returns:
    -----------
        - Data frame in LIONESS format, with rows as edges and columns as samples
    '''
    
    # Create data frames from input files
    pan = pd.read_csv(panda, sep = " ", engine = "python")
    pan["TF-gene"] = "TF_" + pan["tf"] + "<==>" + pan["gene"]

    # Lioness file does not have any header or column names, needs them for t-test later    
    samps = []
    with open(names, 'r') as file:
        for line in file: 
            samps.append(line.strip())
            
    lion.index = pan["TF-gene"] 

    if start is not None:
        lion.columns = samps[start-1:end]        
    else:
        lion.columns = samps 
    
    lion.to_pickle(outfile)   
    
    return(lion)
 
# For hot-starting in case of crashing, so the user does not need to reconstruct all the networks again
# Note that this is just a temporary fix and will be better-implemented in the future.
if __name__ == "__main__":
    parser = argparse.ArgumentParser()    
    parser.add_argument('-p', '--panda', help='Path to panda output file')    
    parser.add_argument('-l', '--lion', help='lioness data frame')    
    parser.add_argument('-t', '--type', help='file type of lioness file, either npy or txt')
    parser.add_argument('-n', '--names', help='File with list of sample names (one per line) in the same order that were supplied for panda/lioness')
    parser.add_argument('-o', '--outfile', help='Path to output file in pickle format (e.g. lioness.pickle)')
    parser.add_argument('-s', '--start', default=None, help='OPTIONAL, the starting index of samples that networks were created for. Only specify this value if "start" and "stop" params were used in the params.yml file.')
    parser.add_argument('-e', '--end', default=None, help='Prints the version of SiSaNA currently being used.')
    args = parser.parse_args()
    
    lioness_array = np.load(args.lion)
    lioness_df = pd.DataFrame(lioness_array)
    
    convert_lion_to_pickle(args.panda, 
                           lioness_df,
                           args.type,
                           args.names,
                           args.outfile,
                           args.start,
                           args.end)    
