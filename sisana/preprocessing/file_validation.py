# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
from pathlib import Path
import sys

def validate_user_params(params_dict, command, subcommand=None):
    """
    Description:
        This code takes the user input params file and validates it for being in the correct format
        to help prevent downstream issues
        
    Parameters:
    -----------
        - params_dict: Dict, the params dictionary read in previously from the user's params file
        - command: the command the user is running (preprocess, generate, etc.)
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
    params["compare"]["required"] = ["datafile", "mapfile", "groups", "filetype"]
    params["compare"]["optional"] = ["datatype", "testtype", "rankby", "outdir"]
    params["compare"]["example"] = """
    compare: 
        datafile: ./example_inputs/lioness_df_indegree_3_decimal_places_subset_200_LumALumB_samps.csv 
        mapfile: ./example_inputs/BRCA_TCGA_200_LumA_LumB_samps_mapping_w_header.csv
        datatype: degree 
        groups:  
            - LumA
            - LumB
        testtype: mw 
        filetype: csv 
        rankby: mediandiff 
        outdir: ./output/compare_means/"""
    
    params["survival"] = {}
    params["survival"]["required"] = ["metadata", "filetype", "sampgroup_colname", "alivestatus_colname", "days_colname", "groups"]
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
    params["gsea"]["optional"] = ["outdir"]
    params["gsea"]["example"] = """
    gsea:
        genefile: ./output/compare_means/comparison_mw_between_LumA_LumB_degree_ranked_mediandiff.rnk
        gmtfile: ./example_inputs/c2.cp.reactome.v2023.2.Hs.symbols.gmt 
        geneset: Reactome 
        outdir: ./output/gsea/"""
    
    params["volcano"] = {}
    params["volcano"]["required"] = ["statsfile", "diffcol"]
    params["volcano"]["optional"] = ["adjpcol", "xaxisthreshold", "adjpvalthreshold", "difftype", "outdir", "top", "numlabels", "genelist"]
    params["volcano"]["example"] = """
    visualize:
        volcano: 
            statsfile: ./output/compare_means/comparison_mw_between_LumA_LumB_degree.txt 
            diffcol: difference_of_medians_(LumB-LumA) 
            adjpcol: FDR 
            xaxisthreshold: 50 
            adjpvalthreshold: 0.25
            difftype: median 
            outdir: ./output/volcano/"""
    
    params["quantity"] = {}
    params["quantity"]["required"] = ["datafile", "statsfile", "filetype", "metadata", "groups", "colors", "prefix", "yaxisname"]    
    params["quantity"]["optional"] = ["plottype", "outdir", "prefix", "genelist", "numgenes", "top"]    
    params["quantity"]["example"] = """
    visualize:
        quantity: 
            datafile: ./example_inputs/lioness_df_indegree_3_decimal_places_subset_200_LumALumB_samps.csv 
            statsfile: ./output/compare_means/comparison_mw_between_LumA_LumB_degree.txt 
            filetype: csv 
            metadata: ./example_inputs/BRCA_TCGA_200_LumA_LumB_samps_mapping_w_header.csv 
            plottype: boxplot 
            groups: 
                - LumA
                - LumB    
            colors: 
                - cornflowerblue
                - orange        
            prefix: LumA_LumB_indegree 
            yaxisname: Indegree
            outdir: ./output/plot_quantity/"""
    
    params["heatmap"] = {}
    params["heatmap"]["required"] = ["datafile", "filetype", "statsfile", "metadata", "genelist", "category_label_columns", "category_column_colors"]
    params["heatmap"]["optional"] = ["column_cluster", "row_cluster", "plot_gene_names", "plot_sample_names", "outdir", "prefix"]
    params["heatmap"]["example"] = """
    visualize:
        heatmap: 
            datafile: ./example_inputs/lioness_df_indegree_3_decimal_places_subset_200_LumALumB_samps.csv 
            filetype: csv 
            statsfile: ./output/compare_means/comparison_mw_between_LumA_LumB_degree.txt 
            metadata: ./example_inputs/BRCA_TCGA_200_LumA_LumB_samps_mapping_w_header.csv 
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
    params["extract"]["required"] = ["pickle", "sampnames", "outdir"]    
    params["extract"]["example"] = """
    extract:
        pickle: ./tmp/lioness.pickle
        sampnames: ./tmp/samples.txt 
        symbols: ./example_inputs/genes_to_extract.txt
        outdir: ./output/extract/"""

    # Ensure the commands are in the params file. If using a visualize command, then you also need 
    # to specify a subcommand in the params file
    if command != "visualize":
        if command not in list(params_dict.keys()):
            raise Exception(f"Error: There is no {command} command present in your params file. Please ensure it is there and that the spelling is correct.")
    else:
        if "visualize" not in list(params_dict.keys()):
            raise Exception(f"Error: There is no visualize command present in your params file. Please ensure it is there and that the spelling is correct.")
        if subcommand not in list(params_dict["visualize"].keys()):
            raise Exception(f"Error: There is no {subcommand} command present in your params file. Please ensure it is there and that the spelling is correct.")
        
    def _validate_required_params(user_params, com, required_params_list, optional_params_list, subcommand=None):
        """
        Description:
            Takes a command name and checks to make sure the required parameters were supplied by the user for that command, as well as ensures that the optional paramaters supplied are valid names
        
        - user_params: dict, a dictionary of the parameters the user has set in their yaml file
        - com: str, the subcommand the user is running (preprocess, generate, etc.)
        - required_params_list: list, list of required parameters for that command
        - optional_params_list: list, list of optional parameters for that command
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
    
    Please ensure they there and their spelling is correct.
    
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
    if command == "generate":
        _validate_required_params(params_dict, "generate", params["generate"]["required"], params["generate"]["optional"])
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
            
    print("Params file structure appears valid. Continuing...")
    
def validate_header(df, delim):
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
            print("number found")
            return True
        except ValueError:
            return False
        
    for i in header_list:
        if _is_number(i):    
            raise Exception(f"Your expression file has numeric values in the column names. Please ensure only strings are supplied in the header of the file. The value in the header that was found which has caused this error is {i}")
    
    print(f"Header of data file appears valid. Continuing...")
           
def validate_metadata(df):
    mapfile = pd.read_csv(df, index_col=0)
    
    if len(mapfile.columns) > 2:
        raise Exception("Error: Please only supply two columns for your mapping file. Ensure the first column is the name and the second column is the group. Also ensure that the file has a header.")

    unique_groups = mapfile.iloc[:, 0].unique()
    if mapfile.columns[0] in unique_groups:
        raise Exception("Error: It appears you do not have a header in your mapping file. Please supply a header with contents that are unique from the values in the columns.")
    
    print(f"Header of metadata file appears valid. Continuing...")
    
def check_genelist_top(user_params, updated_params_w_def, com):
    """
    Description:
        Checks to make sure that at least one of the two options are set for the visualize commands that can either plot
        just the top genes or take a list of genes to plot
     
    - user_params: dict, the params dict from the user's params yaml file 
    - updated_params_w_def: dict, the updated params dict (after filling in non-required default values)   
    - com: str, the subcommand the user is running (preprocess, generate, etc.)
    """
    
    # print(user_params["visualize"][com])
    # print("\n")
    # print(updated_params_w_def["visualize"][com])
    try:
        genelist_user_value = user_params["visualize"][com]["genelist"]
    except KeyError:
        genelist_user_value = None   
         
    try:
        top_user_value = user_params["visualize"][com]["top"]
        top_user_value_supplied = True
    except KeyError:
        top_user_value = False
        top_user_value_supplied = False
   
    # genelist_updated_value = updated_params_w_def["visualize"][com]["genelist"]
    # top_updated_value = updated_params_w_def["visualize"][com]["top"]

    # print("\nOriginal")
    # print(genelist_user_value)
    # print(top_user_value)
    
    # print("\nUpdated")
    # print(genelist_updated_value)
    # print(top_updated_value)
    
    # If user did supply a value for genelist, but then did not submit one for top, then we need to account for the 
    # fact that the code will automatically update the missing 'top' value to True. Doing so would cause issues in the 
    # if statements below, since it would interpret it as the user trying to incorrectly input both a genelist value
    # and use only the top genes
    # 
    # Note that the default value for top is True in the default_parameters.py script, so in this case a user who did not
    # supply a top parameter would automatically have the top_updated_value set to True
        
    if genelist_user_value is not None and top_user_value_supplied == True and top_user_value == True:    
        raise Exception("Error: You have set a value for your 'genelist' and also tried to plot the top values by setting 'top' to True. Only one can be used at a time.")
        
    # if genelist_user_value_supplied is not None and top_user_value_supplied == False:
    # top_updated_value != top_user_value:
    #     elif genelist_user_value is not None and top_user_value == True:
        
    if genelist_user_value == None and top_user_value == False:
        raise Exception("Error: Please make sure to set values for either the 'genelist' parameter or 'top' parameter in your params yaml file.")
      
    # print(genelist_user_value)
    # print(top_user_value)
    
    print("Visualization parameters appear fine. Continuing...")
    # sys.exit(0)
