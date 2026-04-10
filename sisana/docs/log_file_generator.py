import sys
import os

def create_log_file(subcommand: str, params_dict: dict, filenames: list, netzoopy_version: str, sisana_version: str, additional_info=None) -> None: 
    """
    Description:
        This function creates log files for each command the user performs
        
    Parameters:
    -----------
        - subcommand: str, The name of the subcommand used
        - params_dict: dict, The dictionary of the parameters the user supplied
        - filenames: list, The list of file paths that were created with the subcommand
        - additional_info: Dictionary of additional key-value pairs to add to the end of the log file
        
    Returns:
    -----------
        - Nothing
    """
    os.makedirs("./log_files/", exist_ok=True)
    
    basename = f"{subcommand}_log.txt"
    file_outloc = os.path.join("./log_files/", basename)
    
    # Remove the "./" prefix to file names if found. This is for sake of 
    # clarity, since otherwise some file names had them and some did not,
    # just depending on how the user defined them in params file
    fixed_names = [n[2:] if n[:2] == "./" else n for n in filenames]
    
    with open(file_outloc, "w") as file:
        file.write(f"Analysis directory: {os.getcwd()} \n")
        file.write(f"SiSaNA version: {sisana_version} \n")
        file.write(f"netZooPy version: {netzoopy_version} \n")
        
        file.write("\nParameters used:\n")
        for k,v in params_dict.items():
            file.write(f"  {k}: {v}\n")
        
        file.write("\nFiles generated:\n")
        for i in fixed_names: 
            file.write("  - " + i + "\n")

        file.write("\nAdditional information:\n") 
        if additional_info is not None:      
            for k,v in additional_info.items():
                file.write(f"  {k}: {v}\n")

            
