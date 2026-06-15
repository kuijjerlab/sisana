# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
from pathlib import Path
import sys
from sisana.exceptions import NumColorsNumGroupsMismatchError, TooManyCoresError, ExtensionMismatchError

def validate_user_params(params_dict, command, subcommand=None):
    """
    Description:
        This code takes the user input params file and validates it for being in the correct format
        to help prevent downstream issues
        
    Parameters:
    -----------
        - params_dict: Dict, the params dictionary read in previously from the user's params file
        - command: the command the user is running (preprocess, reconstruct, etc.)
        - subcommand: the command the user is running (only if running the visualize command)
    Returns:
    -----------
        - Nothing
    """
    
    params = {}
    
    params["preprocess"] = {}
    params["preprocess"]["required"] = ["exp_file", "filetype"]
    params["preprocess"]["optional"] = ["number", "outdir"]
    params["preprocess"]["example"] = """
    preprocess:
        exp_file: ./example_inputs/BRCA_TCGA_20_LumA_LumB_samps_5000_genes_exp.tsv 
        filetype: tsv 
        number: 5 
        outdir: ./output/preprocess"""    
    
    params["reconstruct"] = {} 
    params["reconstruct"]["required"] = ["exp", "motif", "ppi"]
    params["reconstruct"]["optional"] = ["method", "modeProcess", "pandafilepath", "compute", "ncores", "lionessfilepath", "start", "end"]
    params["reconstruct"]["example"] = """
    reconstruct:
        exp: ./output/preprocess/BRCA_TCGA_20_LumA_LumB_samps_5000_genes_exp_preprocessed.txt
        motif: ./example_inputs/motif_prior_names_2024.tsv
        ppi: ./example_inputs/ppi_prior_2024.tsv 
        method: lioness
        modeProcess: intersection 
        pandafilepath: ./output/network/panda_network.txt 
        compute: cpu 
        ncores: 20
        lionessfilepath: ./output/network/lioness_networks.npy"""

    params["generate"] = {} 
    params["generate"]["required"] = ["exp", "motif", "ppi"]
    params["generate"]["optional"] = ["method", "modeProcess", "pandafilepath", "compute", "ncores", "lionessfilepath", "start", "end"]
    params["generate"]["example"] = """
    generate:
        exp: ./output/preprocess/BRCA_TCGA_20_LumA_LumB_samps_5000_genes_exp_preprocessed.txt
        motif: ./example_inputs/motif_prior_names_2024.tsv
        ppi: ./example_inputs/ppi_prior_2024.tsv 
        method: lioness
        modeProcess: intersection 
        pandafilepath: ./output/network/panda_network.txt 
        compute: cpu 
        ncores: 20
        lionessfilepath: ./output/network/lioness_networks.npy"""
    
    params["compare"] = {}
    params["compare"]["required"] = ["datafile", "mapfile", "groups", "data_filetype", "map_filetype"]
    params["compare"]["optional"] = ["datatype", "testtype", "rankby", "outdir"]
    params["compare"]["example"] = """
    compare: 
        datafile: ./example_inputs/lioness_df_indegree_3_decimal_places_subset_200_LumALumB_samps.csv 
        data_filetype: csv 
        mapfile: ./example_inputs/BRCA_TCGA_200_LumA_LumB_samps_mapping_w_header.csv
        map_filetype: csv
        datatype: degree 
        groups:  
            - LumA
            - LumB
        testtype: mw 
        rankby: mediandiff 
        outdir: ./output/compare_means/"""
    
    params["survival"] = {}
    params["survival"]["required"] = ["metadata", "filetype", "sampgroup_colname", "alivestatus_colname", "days_colname", "groups", "colors"]
    params["survival"]["optional"] = ["outdir"]
    params["survival"]["example"] = """
    survival:
        metadata: ./example_inputs/BRCA_TCGA_200_LumA_LumB_samps_survival_data.csv
        filetype: csv 
        sampgroup_colname: PAM50_subtype 
        alivestatus_colname: Survival_status(False_equals_alive)
        days_colname: days_to_death_or_last_followup 
        groups:
            - LumA
            - LumB
        outdir: ./output/survival/"""
    
    params["gsea"] = {}
    params["gsea"]["required"] = ["genefile", "gmtfile", "geneset"]
    params["gsea"]["optional"] = ["color", "outdir"]
    params["gsea"]["example"] = """
    gsea:
        genefile: ./output/compare_means/comparison_mw_between_LumA_LumB_degree_ranked_mediandiff.rnk
        gmtfile: ./example_inputs/c2.cp.reactome.v2023.2.Hs.symbols.gmt 
        geneset: Reactome 
        outdir: ./output/gsea/"""
    
    params["volcano"] = {}
    params["volcano"]["required"] = ["statsfile", "diffcol", "top", "genelist", "groups", "colors"]
    params["volcano"]["optional"] = ["adjpcol", "xaxisthreshold", "adjpvalthreshold", "difftype", "outdir", "numlabels"]
    params["volcano"]["example"] = """
    visualize:
        volcano: 
            statsfile: ./output/compare_means/comparison_mw_between_LumA_LumB_degree.txt 
            diffcol: difference_of_medians_(LumB-LumA) 
            adjpcol: FDR 
            xaxisthreshold: 50 
            adjpvalthreshold: 0.25
            difftype: median 
            outdir: ./output/volcano/
            genelist: ./example_inputs/volcano_plot_genes.txt
            top: False"""
    
    params["quantity"] = {}
    params["quantity"]["required"] = ["datafile", "statsfile", "data_filetype", "metadata", "groups", "colors", "prefix", "yaxisname",  "genelist", "top", "metadata_filetype"]    
    params["quantity"]["optional"] = ["plottype", "outdir", "prefix", "numgenes"]    
    params["quantity"]["example"] = """
    visualize:
        quantity: 
            datafile: ./example_inputs/lioness_df_indegree_3_decimal_places_subset_200_LumALumB_samps.csv 
            statsfile: ./output/compare_means/comparison_mw_between_LumA_LumB_degree.txt 
            data_filetype: csv 
            metadata: ./example_inputs/BRCA_TCGA_200_LumA_LumB_samps_mapping_w_header.csv
            metadata_filetype: csv 
            plottype: boxplot 
            groups: 
                - LumA
                - LumB    
            colors: 
                - cornflowerblue
                - orange        
            prefix: LumA_LumB_indegree 
            yaxisname: Indegree
            outdir: ./output/plot_quantity/
            genelist: ./example_inputs/quantity_plot_genes.txt
            top: False"""
    
    params["heatmap"] = {}
    params["heatmap"]["required"] = ["datafile", "data_filetype", "statsfile", "metadata", "genelist", "category_label_columns", "category_column_colors", "metadata_filetype"]
    params["heatmap"]["optional"] = ["column_cluster", "row_cluster", "plot_gene_names", "plot_sample_names", "outdir", "prefix", "subset_for", "data_color"]
    params["heatmap"]["example"] = """
    visualize:
        heatmap: 
            datafile: ./example_inputs/lioness_df_indegree_3_decimal_places_subset_200_LumALumB_samps.csv 
            data_filetype: csv 
            statsfile: ./output/compare_means/comparison_mw_between_LumA_LumB_degree.txt 
            metadata: ./example_inputs/BRCA_TCGA_200_LumA_LumB_samps_mapping_w_header.csv 
            metadata_filetype: csv
            genelist: ./example_inputs/heatmap_genes.txt 
            column_cluster: False
            row_cluster: True 
            plot_gene_names: True 
            plot_sample_names: False 
            category_label_columns:
                - group
            category_column_colors: 
                - {'LumA': 'cornflowerblue', 
                  'LumB': 'orange'}
            outdir: ./output/heatmap/
            prefix: TCGA_200_LumA_LumB_samps """
    
    params["extract"] = {}
    params["extract"]["required"] = ["symbols"]    
    params["extract"]["optional"] = ["pickle", "sampnames", "outdir"]    
    params["extract"]["example"] = """
    extract:
        pickle: ./tmp/lioness.pickle
        sampnames: ./tmp/samples.txt 
        symbols: ./example_inputs/genes_to_extract.txt
        outdir: ./output/extract/"""    
    
    params["quantnorm"] = {}
    params["quantnorm"]["required"] = ["network_file", "filetype"]    
    params["quantnorm"]["optional"] = ["pandafilepath", "num_cpus", "start", "end"]    
    params["quantnorm"]["example"] = """
    quantnorm:
        network_file: ./tmp/lioness.npy
        filetype: npy
        pandafilepath: ./output/network/panda_network.txt
        num_cpus: 4"""

    # Ensure the commands are in the params file. If using a visualize command, then you also need 
    # to specify a subcommand in the params file
    if command != "visualize":
        if command not in list(params_dict.keys()):
            raise Exception(f"Error: There is no {command} command present in your params file. Please ensure it is there and that the spelling is correct.")
    else:
        if subcommand is None:
            raise Exception(f"Error: You are running the visualize command, but you have not specified a subcommand (heatmap, volcano, or quantity) in your params file. Please ensure you have specified one of these subcommands and that the spelling is correct.")
        if "visualize" not in list(params_dict.keys()):
            raise Exception(f"Error: There is no visualize command present in your params file. Please ensure it is there and that the spelling is correct.")
        if subcommand not in list(params_dict["visualize"].keys()):
            raise Exception(f"Error: There is no {subcommand} command present in your params file. Please ensure it is there and that the spelling is correct.")
        
    def _validate_required_params(user_params, com, required_params_list, optional_params_list, subcommand=None):
        """
        Description:
            Takes a command name and checks to make sure the required parameters were supplied by the user for that command, as well as ensures that the optional paramaters supplied are valid names
        
        Parameters:
        -----------     
            - user_params: dict, a dictionary of the parameters the user has set in their yaml file
            - com: str, the subcommand the user is running (preprocess, reconstruct, etc.)
            - required_params_list: list, list of required parameters for that command
            - optional_params_list: list, list of optional parameters for that command
        
        Returns:
        -----------
            - Nothing           
        """

        if command != "visualize":
            command_dict = user_params[com]
        else:
            command_dict = user_params["visualize"][com]
        
        given_commands = list(command_dict.keys())
        
        # Check if required params are given in the user's param file
        not_supplied_req_params = []
        for i in required_params_list:
            if i not in given_commands:
                not_supplied_req_params.append(i)
                
        if not_supplied_req_params: # if list is not empty
            raise Exception(f"""
                    Error: You are missing the following required {com} parameters in your params file: {', '.join(not_supplied_req_params)}
                    
                    Please ensure they are there and their spelling is correct.
                    
                    Please see the following for an example:
                    {params[com]["example"]}""")
            
        # Check to see if optional params are valid param options
        unrecognized_opt_params = []
        for i in given_commands:
            # print(i)
            if (i not in optional_params_list) and (i not in required_params_list):
                unrecognized_opt_params.append(i)
        
        if unrecognized_opt_params: # if list is not empty
            raise Exception(f"Error: the following optional {com} parameters you supplied in your params file are not recognized: {', '.join(unrecognized_opt_params)}")
    
    if command == "preprocess":
        _validate_required_params(params_dict, "preprocess", params["preprocess"]["required"], params["preprocess"]["optional"])
    if command == "reconstruct" or command == "generate":
        _validate_required_params(params_dict, command, params[command]["required"], params[command]["optional"])
    if command == "compare":
        _validate_required_params(params_dict, "compare", params["compare"]["required"], params["compare"]["optional"])
    if command == "survival":
        _validate_required_params(params_dict, "survival", params["survival"]["required"], params["survival"]["optional"])
    if command == "gsea":
        _validate_required_params(params_dict, "gsea", params["gsea"]["required"], params["gsea"]["optional"])
    if subcommand == "volcano":
        _validate_required_params(params_dict, "volcano", params["volcano"]["required"], params["volcano"]["optional"])
    if subcommand == "quantity":
        _validate_required_params(params_dict, "quantity", params["quantity"]["required"], params["quantity"]["optional"])
    if subcommand == "heatmap":
        _validate_required_params(params_dict, "heatmap", params["heatmap"]["required"], params["heatmap"]["optional"])
        # _check_genelist_top(params_dict, "heatmap") # Commenting out for now since the "top" option is still being implemented for the heatmap
    if command == "extract":
        _validate_required_params(params_dict, "extract", params["extract"]["required"], params["extract"]["optional"])
            
    print("Params file contains validly named parameters. Continuing...")
    
def check_for_header(df, delim):
    """
    Description:
        Checks to make sure that the supplied file has a header

    Parameters:
    -----------     
        - requested_cores: int, the number of cores the user requested in the params.yml file
    
    Returns:
    -----------
        - Nothing
    """
    
    # Make sure header is supplied (e.g. no numeric vals)
    if delim == "csv":
        with open(df, 'r') as f:
            header = f.readline()
            header_list = header.split(",")[1:]
    else:
        with open(df, 'r') as f:
            header = f.readline()
            header_list = header.split()[1:]
    
    # Check if each value in the header is a string. Note that they are automatically saved as a string when read in by the reader previously so we need to attempt to convert
    # to a numeric value. Ints and floats will convert fine, raising the error.
    def _is_number(s):
        try:
            float(s)
            return True
        except ValueError:
            return False
        
    for i in header_list:
        if _is_number(i):    
            raise Exception(f"Your expression file has numeric values in the column names. Please ensure only strings are supplied in the header of the file. The value in the header that was found which has caused this error is {i}")
    
    print(f"Header of data file appears valid. Continuing...")
           
def check_ncore_value(requested_cores):
    """
    Description:
        Checks to make sure that the number of cores requested by the user is not more than the number of samples. This is necessary
        not only as a way to reduce computational load, but also because Lioness will encounter the error "AttributeError: 'Lioness' object 
        has no attribute 'total_lioness_network'" if not enforced.

    Parameters:
    -----------     
        - requested_cores: int, the number of cores the user requested in the params.yml file
    
    Returns:
    -----------
        - Nothing
    """
    
    with open('./tmp/num_samples.txt') as f:
        nsamps = int(f.read())
            
    if int(requested_cores) > nsamps:
        raise TooManyCoresError(requested_cores, nsamps)

def validate_metadata(df: pd.DataFrame, testtype: str, groups: list=None):
    """
    Description:
        Checks the metadata for correct formatting

    Parameters:
    -----------     
        - df: pd.DataFrame, Metadata data frame that the user supplied in the params file
        - testtype: str, the type of test the user is running (tt, mw, paired_tt, wilcoxon)
        - groups: list, the two groups the user is comparing (only required for paired tests)
    
    Returns:
    -----------
        - Nothing
    """

    if testtype == "paired_tt" or testtype == "wilcoxon":
        colnames = [df.index.name, df.columns[0]]
        print(colnames)
        print(groups)
        
        if (sorted(colnames) != sorted(groups)):
            raise Exception(f"Error: You have specified that you are running a paired test, but your specified groups in your mapping file do not match the column names of your metadata file.")
            
    else:
        if df.columns[0] in df.iloc[:, 0].unique():
            raise Exception("Error: It appears you do not have a header in your mapping file. Please supply a header with contents that are unique from the values in the columns.")
    
    print(f"Header of metadata file appears valid. Continuing...")
    
def check_genelist_top(user_params, updated_params_w_def, com):
    """
    Description:
        Checks to make sure that at least one of the two options are set for the visualize commands that can either plot
        just the top genes or take a list of genes to plot
     
    Parameters:
    -----------  
        - user_params: dict, the params dict from the user's params yaml file 
        - updated_params_w_def: dict, the updated params dict (after filling in non-required default values)   
        - com: str, the subcommand the user is running (preprocess, reconstruct, etc.)
        
    Returns:
    -----------
        - Nothing
    """
    if user_params["visualize"][com]["genelist"] is None and user_params["visualize"][com]["top"] == False:
        raise Exception("Error: You must choose to either plot a list of genes by supplying a file to the 'genelist' parameter or plot only the top genes by setting the 'top' parameter to True.")
    
    if user_params["visualize"][com]["genelist"] is not None and user_params["visualize"][com]["top"] == True:
        raise Exception("Error: You must choose to either plot a list of genes by supplying a file to the 'genelist' parameter or plot only the top genes by setting the 'top' parameter to True. Only one can be used at a given time.")
    
    # # print(user_params["visualize"][com])
    # # print("\n")
    # # print(updated_params_w_def["visualize"][com])
    # try:
    #     genelist_user_value = user_params["visualize"][com]["genelist"]
    # except KeyError:
    #     genelist_user_value = None   
         
    # try:
    #     top_user_value = user_params["visualize"][com]["top"]
    #     top_user_value_supplied = True
    # except KeyError:
    #     top_user_value = False
    #     top_user_value_supplied = False
   
    # # genelist_updated_value = updated_params_w_def["visualize"][com]["genelist"]
    # # top_updated_value = updated_params_w_def["visualize"][com]["top"]

    # # print("\nOriginal")
    # # print(genelist_user_value)
    # # print(top_user_value)
    
    # # print("\nUpdated")
    # # print(genelist_updated_value)
    # # print(top_updated_value)
    
    # # If user did supply a value for genelist, but then did not submit one for top, then we need to account for the 
    # # fact that the code will automatically update the missing 'top' value to True. Doing so would cause issues in the 
    # # if statements below, since it would interpret it as the user trying to incorrectly input both a genelist value
    # # and use only the top genes
    # # 
    # # Note that the default value for top is True in the default_parameters.py script, so in this case a user who did not
    # # supply a top parameter would automatically have the top_updated_value set to True
        
    # if genelist_user_value is not None and top_user_value_supplied == True and top_user_value == True:    
    #     raise Exception("Error: You have set a value for your 'genelist' and also tried to plot the top values by setting 'top' to True. Only one can be used at a time.")
        
    # # if genelist_user_value_supplied is not None and top_user_value_supplied == False:
    # # top_updated_value != top_user_value:
    # #     elif genelist_user_value is not None and top_user_value == True:
        
    # if genelist_user_value == None and top_user_value == False:
    #     raise Exception("Error: Please make sure to set values for either the 'genelist' parameter or 'top' parameter in your params yaml file.")
      
    # # print(genelist_user_value)
    # # print(top_user_value)
    
    # print("Visualization parameters appear fine. Continuing...")
    # # sys.exit(0)
    
def check_no_hyphens_in_group_names(meta: str):
    """
    Description:
        Checks to make sure that none of the group names contain hyphens, as this can cause issues with downstream analysis.
     
    Parameters:
    -----------  
        - meta: str, the metadata df
        
    Returns:
    -----------
        - Nothing
    """
    
    for groupname in list(meta.iloc[:, 0]):
        if "-" in groupname:
            raise ValueError(f"Group name '{groupname}' contains a hyphen, which is not allowed.")
    print("Group names in mapping file appear valid. Continuing...")
    
def check_num_group_colors(user_params: dict, com: str):
    """
    Description:
        Checks to make sure that the number of groups in the metadata file matches the number of color codes provided.
     
    Parameters:
    -----------  
        - user_params: dict, the user-defined parameters for the visualization step.
        - com: str, the subcommand the user is running (e.g. volcano, quantity, heatmap)
        
    Returns:
    -----------
        - Nothing
    """
    if len(user_params["visualize"][com]["groups"]) != len(user_params["visualize"][com]["colors"]):
        raise NumColorsNumGroupsMismatchError("groups", "colors", len(user_params["visualize"][com]["groups"]), len(user_params["visualize"][com]["colors"]))
    
def check_file_extension(filename, ext):
    """
    Description:
        Checks whether the user-supplied file has the correct file extension based on the delimiter they specified in the params file. For example, if they specified a comma delimiter, 
        then the file should have a .csv extension. If they specified a tab delimiter, then the file should have a .txt or .tsv extension.

    Parameters:
    -----------     
        - filename: str, the path to the file
        - ext: str, the user-supplied extension 
    
    Returns:
    -----------
        - Nothing
    """
    
    if ext == "csv" and (filename[-3:] == "txt" or filename[-3:] == "tsv"):
        raise ExtensionMismatchError(filename, ext)
    if (ext == "txt" or ext == "tsv") and filename[-3:] == "csv":
        raise ExtensionMismatchError(filename, ext)
