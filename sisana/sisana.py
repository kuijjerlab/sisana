import yaml
import argparse
from importlib.metadata import version
from sisana.default_parameters import get_default_params 
from sisana.exceptions import *
from sisana.preprocessing import preprocess_data, validate_user_params, check_for_header, validate_metadata, check_genelist_top, check_ncore_value, check_no_hyphens_in_group_names, check_num_group_colors
from sisana.postprocessing import convert_lion_to_pickle, extract_tfs_genes, combine_files, quantile_normalize_edges
from sisana.analyze_networks import calculate_panda_degree, calculate_lioness_degree, compare_bw_groups, survival_analysis, perform_gsea, plot_volcano, plot_expression_degree, plot_heatmap, plot_clustermap, summarize
from sisana.reconstruct import make_panda_network, make_lioness_networks
from sisana.example_input import find_example_paths, fetch_files
import sisana.docs
from sisana.docs import create_log_file
import os 
import pandas as pd
import sys
import re
import glob
import numpy as np
from pathlib import Path
 
def cli():
    """
    SiSaNA command line interface
    """

    DESCRIPTION = """
    SiSaNA - Single Sample Network Analysis
    A command line interface tool used to reconstruct and analyze 
    PANDA and LIONESS networks. It works through subcommands. 
    The command 'sisana reconstruct params.yaml', for example,
    will reconstruct a PANDA or LIONESS network, using the parameters 
    set in the params.yaml file.
    Developed by Nolan Newman (nolan.newman@ncmm.uio.no).
    """
    EPILOG = """
    Code available under MIT license:
    https://github.com/kuijjerlab/sisana
    """
    
    # create the top-level parser
    parser = argparse.ArgumentParser(prog='sisana.py', description=DESCRIPTION, epilog=EPILOG)    
    parser.add_argument('-e', '--example', action='store_true', help='Flag; Copies the example input files into a directory called "./example_inputs"')    
    parser.add_argument('-s', '--setAndForget', action='store_true', help='Flag; Will attempt to run ALL STEPS of SiSaNA at once. Warning: This requires a very well-formatted params file and should not be used by first-time users. Most users will want to run each of the steps individually."')    
    parser.add_argument('-v', '--version', action='store_true', help='Prints the version of SiSaNA currently being used.')

    # Add subcommands
    subparsers = parser.add_subparsers(title='Subcommands', dest='command')
    pre = subparsers.add_parser('preprocess', help='Filters expression data for parameters (e.g. genes) that are only present in at least m samples. Also filters each input file so they have the same genes and TFs across each', epilog=sisana.docs.preprocess_desc, formatter_class=argparse.RawDescriptionHelpFormatter)
    rec = subparsers.add_parser('reconstruct', aliases=['generate'], help='reconstructs PANDA and LIONESS networks', epilog=sisana.docs.reconstruct_desc, formatter_class=argparse.RawDescriptionHelpFormatter)
    comb = subparsers.add_parser('combine', help='Combines indegree and outdegree files ran in batches', epilog=sisana.docs.combine_desc, formatter_class=argparse.RawDescriptionHelpFormatter)
    ext = subparsers.add_parser('extract', help='Extract edges connected to specified TFs/genes', epilog=sisana.docs.extract_desc, formatter_class=argparse.RawDescriptionHelpFormatter)
    comp = subparsers.add_parser('compare', help='Compare networks between sample groups', epilog=sisana.docs.compare_desc, formatter_class=argparse.RawDescriptionHelpFormatter)
    surv = subparsers.add_parser('survival', help='Compare survival times of individuals between sample groups', epilog=sisana.docs.survival_desc, formatter_class=argparse.RawDescriptionHelpFormatter)
    gsea = subparsers.add_parser('gsea', help='Perform gene set enrichment analysis between sample groups', epilog=sisana.docs.gsea_desc, formatter_class=argparse.RawDescriptionHelpFormatter)
    vis = subparsers.add_parser('visualize', help='Visualize the calculated degrees of each sample group', epilog=sisana.docs.visualize_desc, formatter_class=argparse.RawDescriptionHelpFormatter)
    summ = subparsers.add_parser('summarize', aliases=["summarise"], help='Summarize the outputs in an html file', epilog=sisana.docs.summarize_desc, formatter_class=argparse.RawDescriptionHelpFormatter)
    qnorm = subparsers.add_parser('quantnorm', help='Quantile normalize network edges', epilog=sisana.docs.quantnorm_desc, formatter_class=argparse.RawDescriptionHelpFormatter)

    # options for preprocess subcommand
    pre.add_argument("params", type=str, help='Path to yaml file containing the parameters to use')
        
    # options for reconstruct subcommand    
    rec.add_argument("params", type=str, help='Path to yaml file containing the parameters to use')

    # options for combine subcommand    
    comb.add_argument("params", type=str, help='Path to yaml file containing the parameters to use')

    # options for extract subcommand
    ext.add_argument("extractchoice", type=str, choices = ["genes", "tfs"], help="Do you want to extract specific gene or TF edges?")   
    ext.add_argument("params", type=str, help='Path to yaml file containing the parameters to use')

    # options for compare subcommand
    comp.add_argument("params", type=str, help='Path to yaml file containing the parameters to use')

    # options for survival subcommand
    surv.add_argument("params", type=str, help='Path to yaml file containing the parameters to use')

    # options for gsea subcommand    
    gsea.add_argument("params", type=str, help='Path to yaml file containing the parameters to use')

    # options for visualize subcommand
    vis.add_argument("plotchoice", type=str, choices = ["quantity", "heatmap", "volcano"], nargs='?', help="The type of plot to create")   
    vis.add_argument("params", type=str, help='Path to yaml file containing the parameters to use')
    
    # options for summarize subcommand    
    summ.add_argument("logdir", nargs='?', type=str, default="./log_files/", help='Path to the directory containing the previously made log files')

    # options for quantnorm subcommand
    qnorm.add_argument("params", type=str, help='Path to yaml file containing the parameters to use')

    args = parser.parse_args()

    s_version = sisana.__version__
    nzp_version = version('netZooPy')     
    
    # If user wants example files, retrieve them from Zenodo
    if args.example:
        print("Downloading example input files from Zenodo. Please wait...")
        fetch_files()
        print("Example input files have been created in ./example_inputs/")
        sys.exit(0)
        
    # If user wants version info
    if args.version:
        print(f"SiSaNA version: {s_version}")
        sys.exit(0)

    # If user has already performed analysis and wants an HTML summary file
    if args.command != "summarize" and args.command != "summarise": 
        params = yaml.load(open(args.params), Loader=yaml.FullLoader)
    else:
        summarize(args.logdir)
        sys.exit(0)

    # Validate params file
    if args.command != "visualize":
        validate_user_params(params, args.command)
    else:
        validate_user_params(params, args.command, args.plotchoice)

    # Create output for temp files if one does not already exist
    os.makedirs('./tmp/', exist_ok=True)
    
    # Create a dictionary with the default parameters for each step
    def_params = get_default_params()
                
    def _update_if_different(default_dict, user_dict) -> dict:
        """    
        Description:
            Updates the default_dict with the user-defined parameters from user_dict,
            then returns the resulting default_dict
                    
        Parameters:
        -----------
            - default_dict: dict, a dictionary containing the default sisana parameters
            - user_dict: dict, a dictionary containing the parameters the user defined
            
        Returns:
        -----------
            -  The default dict, with the default parameters updated if the user has supplied
               values for those parameters, otherwise the defaults are kept
        """
        
        temp_dict = default_dict
            
        for key, source_value in user_dict.items():
            if key not in default_dict:
                temp_dict[key] = user_dict[key]
                
            if default_dict[key] != user_dict[key]:
                temp_dict[key] = source_value
    
        return temp_dict

    updated_params = {}
    
    single_dict_keys = ["preprocess", "reconstruct", "generate", "combine", "compare", "survival", "gsea", "extract", "quantnorm"]
    nested_dict_keys = ["volcano", "quantity", "heatmap"]

    for key in single_dict_keys:
        if key in params:
            updated_params[key] = _update_if_different(def_params[key], params[key])

    updated_params["visualize"] = {}
    if "visualize" in params:
        for vis_type in nested_dict_keys:
            if vis_type in params["visualize"]:
                updated_params["visualize"][vis_type] = _update_if_different(def_params["visualize"][vis_type], params["visualize"][vis_type])

    ########################################################
    # 1) Preprocess the data
    ########################################################
    
    if args.command == 'preprocess':
        
        preprocess_params = updated_params['preprocess']
        
        # # Save the order of the sample names to their own file, then export the data frame without a header, since that is what is required for CLI version of PANDA
        # expdf = pd.read_csv(preprocess_params['exp_file'], sep='\t', index_col=0)
        # name_list = list(expdf.columns.values)
        
        # with open('./tmp/samples.txt', 'w') as f:
        #     for line in name_list:
        #         f.write(f"{line}\n")
        
        check_for_header(preprocess_params['exp_file'], preprocess_params['filetype'])
        
        # Remove genes that are not expressed in at least the user-defined minimum ("number")
        results = preprocess_data(preprocess_params['exp_file'], 
                        preprocess_params['filetype'], 
                        preprocess_params['number'],
                        preprocess_params['outdir'])  
        
        fname, genes_kept, genes_removed = results[0], results[1], results[2] 
        
        extra_info_preprocess = {"genes removed": genes_removed, "genes kept": genes_kept}
        
        create_log_file(subcommand="preprocess", 
                        params_dict=preprocess_params, 
                        filenames=[fname], 
                        netzoopy_version=nzp_version,
                        sisana_version=s_version,
                        additional_info=extra_info_preprocess)
                    
    ########################################################
    # 2) Run PANDA/LIONESS, using the parameters from the yaml file
    ########################################################

    if args.command == 'reconstruct' or args.command == 'generate':
        
        reconstruct_params = updated_params[args.command]
        if reconstruct_params["method"] == "lioness":
            check_ncore_value(reconstruct_params["ncores"])
            
        panda_output_location = reconstruct_params["pandafilepath"]
        
        print("Now reconstructing PANDA network...")        
        pan = make_panda_network(reconstruct_params["exp"],
                                    reconstruct_params["motif"],
                                    reconstruct_params["ppi"],
                                    reconstruct_params["compute"],
                                    reconstruct_params["modeProcess"],
                                    pandafilepath=panda_output_location)
            
        # If user wants to run lioness, then we need to do the following
        if reconstruct_params['method'].lower() == 'lioness':
            lion_files = make_lioness_networks(panda=pan,
                                  compute=reconstruct_params["compute"],
                                  ncores=reconstruct_params["ncores"],
                                  start=reconstruct_params["start"],
                                  end=reconstruct_params["end"],
                                  lioness_fpath=reconstruct_params["lionessfilepath"],
                                  panda_fpath=panda_output_location)
            
        print(f"\nPANDA network saved to {panda_output_location}")
        print(f"PANDA degrees saved to:") 
        print(f"{str(panda_output_location)[:-4]}_outdegree.csv")
        print(f"{str(panda_output_location)[:-4]}_indegree.csv")
        
        outfiles = [panda_output_location,
                    f"{str(panda_output_location)[:-4]}_outdegree.csv",
                    f"{str(panda_output_location)[:-4]}_indegree.csv",
                    str(lion_files["lioness_nw_filepath"]),
                    str(lion_files["lioness_indeg_filepath"]),
                    str(lion_files["lioness_outdeg_filepath"])]
            
        create_log_file(subcommand=args.command, 
                        params_dict=reconstruct_params, 
                        netzoopy_version=nzp_version,
                        sisana_version=s_version,
                        filenames=outfiles)
        
    ########################################################
    # 2.5) (OPTIONAL) Combine the multiple degree files into a single
    #      output file. Only used samples were "batched" when creating
    #      the single-sample networks in the previous step
    ########################################################
    if args.command == 'combine':
        
        combine_params = updated_params['combine']
        outfiles = combine_files(combine_params)
        
        create_log_file(subcommand="combine", 
                        params_dict=combine_params, 
                        netzoopy_version=nzp_version,
                        sisana_version=s_version,
                        filenames=outfiles)       

    ########################################################
    # 3) Compare degree (or expression) between sample groups
    ########################################################
        
    if args.command == "compare":     
        compare_means_params = updated_params['compare']
        
        # check_no_hyphens_in_group_names(compare_means_params["mapfile"])        
        # validate_metadata(compare_means_params['mapfile'])

        outfiles = compare_bw_groups(datafile=compare_means_params["datafile"], 
                                    mapfile=compare_means_params["mapfile"], 
                                    datatype=compare_means_params["datatype"], 
                                    groups=compare_means_params["groups"],
                                    testtype=compare_means_params["testtype"], 
                                    data_filetype=compare_means_params["data_filetype"],
                                    map_filetype=compare_means_params["map_filetype"],
                                    rankby_col=compare_means_params["rankby"],
                                    outdir=compare_means_params["outdir"])
 
        create_log_file(subcommand="compare_means", 
                        params_dict=compare_means_params, 
                        netzoopy_version=nzp_version,
                        sisana_version=s_version,
                        filenames=outfiles)
    
    ########################################################
    # 4) Perform gene set enrichment analysis
    ########################################################   
        
    if args.command == 'gsea':    
        gsea_params = updated_params["gsea"]

        outfiles = perform_gsea(genefile=gsea_params["genefile"], 
                        gmtfile=gsea_params["gmtfile"], 
                        geneset=gsea_params["geneset"], 
                        color=gsea_params["color"],
                        outdir=gsea_params["outdir"])
        
        create_log_file(subcommand="gsea", 
                        params_dict=gsea_params, 
                        netzoopy_version=nzp_version,
                        sisana_version=s_version,
                        filenames=outfiles)
                
    ########################################################
    # 5) Visualize results
    ########################################################       

    if args.command == "visualize":                  

        if args.plotchoice == "volcano": 
            # check_genelist_top(params, updated_params, "volcano")
            check_num_group_colors(params, "volcano")
            
            volcano_params = updated_params["visualize"]["volcano"]

            outfiles, down_group, down_gene_count, up_group, up_gene_count = plot_volcano(statsfile=volcano_params["statsfile"],
                         diffcol=volcano_params["diffcol"],
                         adjpcol=volcano_params["adjpcol"],
                         adjpvalthreshold=volcano_params["adjpvalthreshold"],
                         xaxisthreshold=volcano_params["xaxisthreshold"],
                         groups=volcano_params["groups"],
                         colors=volcano_params["colors"],
                         difftype=volcano_params["difftype"],
                         genelist=volcano_params["genelist"],
                         outdir=volcano_params["outdir"],
                         numlabels=volcano_params["numlabels"],
                         top=volcano_params["top"])      
            
            extra_info_num_genes = {"down_group": down_group, "down_gene_count": down_gene_count, "up_group": up_group, "up_gene_count": up_gene_count}
            
            create_log_file(subcommand="volcano_plot", 
                            params_dict=volcano_params, 
                            netzoopy_version=nzp_version,
                            sisana_version=s_version,
                            filenames=outfiles, 
                            additional_info=extra_info_num_genes)
                
        if args.plotchoice == "quantity":  

            check_genelist_top(params, updated_params, "quantity")
 
            try:   
                quantity_params = updated_params["visualize"]["quantity"]
            except KeyError:
                raise Exception("Error: No parameters for visualization of 'quantity' have been set in the params.yml file.")
            
            outfiles = plot_expression_degree(datafile=quantity_params["datafile"],
                        data_filetype=quantity_params["data_filetype"], 
                        statsfile=quantity_params["statsfile"], 
                        metadata=quantity_params["metadata"],
                        metadata_filetype=quantity_params["metadata_filetype"], 
                        plottype=quantity_params["plottype"],
                        groups=quantity_params["groups"],
                        colors=quantity_params["colors"],
                        prefix=quantity_params["prefix"],
                        yaxisname=quantity_params["yaxisname"],
                        outdir=quantity_params["outdir"],
                        genelist=quantity_params["genelist"],
                        top=quantity_params["top"])   
                
            create_log_file(subcommand="quantity_plot", 
                            params_dict=quantity_params, 
                            netzoopy_version=nzp_version,
                            sisana_version=s_version,
                            filenames=outfiles)               
                
        # For now, the plot_heatmap option is being deprecated for use of the plot_clustermap option instead,
        # as the clustermap option allows for more user control and clustering of patients/parameters
        # if args.plotchoice == "heatmap":    
        #     plot_heatmap(datafile=params["visualize"]["heatmap"]["datafile"],
        #                 filetype=params["visualize"]["heatmap"]["filetype"], 
        #                 statsfile=params["visualize"]["heatmap"]["statsfile"],
        #                 metadata=params["visualize"]["heatmap"]["metadata"],
        #                 genelist=params["visualize"]["heatmap"]["genelist"],
        #                 groups=params["visualize"]["heatmap"]["groups"],
        #                 prefix=params["visualize"]["heatmap"]["prefix"],
        #                 plotnames=params["visualize"]["heatmap"]["plotnames"],
        #                 outdir=params["visualize"]["heatmap"]["outdir"],
        #                 top=False)  
            
        if args.plotchoice == "heatmap":  
            # check_genelist_top(params, updated_params, "heatmap")
  
            try:
                heatmap_params = updated_params["visualize"]["heatmap"]
            except KeyError:
                raise Exception("Error: No parameters for visualization of 'heatmap' have been set in the params.yml file.")
            
            outfiles = plot_clustermap(datafile=heatmap_params["datafile"],
                        data_filetype=heatmap_params["data_filetype"], 
                        metadata=heatmap_params["metadata"],
                        metadata_filetype=heatmap_params["metadata_filetype"],
                        genelist=heatmap_params["genelist"],
                        column_cluster=heatmap_params["column_cluster"],
                        row_cluster=heatmap_params["row_cluster"],
                        data_color=heatmap_params["data_color"],
                        prefix=heatmap_params["prefix"],
                        outdir=heatmap_params["outdir"],
                        plot_gene_names=heatmap_params["plot_gene_names"],
                        plot_sample_names=heatmap_params["plot_sample_names"],
                        category_label_columns=heatmap_params["category_label_columns"],
                        category_column_colors=heatmap_params["category_column_colors"],                       
                        top=False,
                        subset_for=heatmap_params["subset_for"])   
            
            create_log_file(subcommand="heatmap", 
                            params_dict=heatmap_params, 
                            netzoopy_version=nzp_version,
                            sisana_version=s_version,
                            filenames=outfiles)  
            
    ########################################################
    # (Optional) Extract edges that connect to specific TFs/genes
    ########################################################

    if args.command == 'extract':
        try:
            extract_params = updated_params["extract"]
        except KeyError:
            raise Exception("Error: No parameters for 'extract' have been set in the params.yml file.")
            
        outfiles = extract_tfs_genes(pickle=extract_params["pickle"], 
                         datatype=args.extractchoice, 
                         sampnames=extract_params["sampnames"],
                         symbols=extract_params["symbols"], 
                         outdir=extract_params["outdir"])
        
        create_log_file(subcommand="extract", 
                        params_dict=extract_params,
                        netzoopy_version=nzp_version, 
                        sisana_version=s_version,
                        filenames=[outfiles])  
                
    ########################################################
    # (Optional) Perform survival analysis
    ########################################################
   
    if args.command == "survival":     
        compare_survival_params = updated_params['survival']

        try:
            outfiles = survival_analysis(metadata=compare_survival_params["metadata"],
                            filetype=compare_survival_params["filetype"], 
                            sampgroup_colname=compare_survival_params["sampgroup_colname"],
                            alivestatus_colname=compare_survival_params["alivestatus_colname"],
                            days_colname=compare_survival_params["days_colname"],
                            groups=compare_survival_params["groups"],
                            colors=compare_survival_params["colors"],
                            outdir=compare_survival_params["outdir"],
                            appendname=compare_survival_params["appendname"])
        except:
            outfiles = survival_analysis(metadata=compare_survival_params["metadata"],
                            filetype=compare_survival_params["filetype"], 
                            sampgroup_colname=compare_survival_params["sampgroup_colname"],
                            alivestatus_colname=compare_survival_params["alivestatus_colname"],
                            days_colname=compare_survival_params["days_colname"],
                            groups=compare_survival_params["groups"],
                            colors=compare_survival_params["colors"],
                            outdir=compare_survival_params["outdir"])
        fnames, pval, sig = outfiles[0], outfiles[1], outfiles[2] 
        
        extra_info = {"p-value": pval, "significant?": sig}
        
        create_log_file(subcommand="compare_survival", 
                        sisana_version=s_version,
                        netzoopy_version=nzp_version,
                        params_dict=compare_survival_params, 
                        filenames=fnames, 
                        additional_info=extra_info)

    # ########################################################
    # # (Optional) TODO: Quantile normalize edges, then calculate degree
    # ########################################################
    if args.command == "quantnorm":  
        qnorm_params = updated_params['quantnorm']

        print(qnorm_params)

        outfiles = quantile_normalize_edges(net=qnorm_params["network_file"], 
                        filetype=qnorm_params["filetype"],
                        pandafilepath=qnorm_params["pandafilepath"],
                        num_cpus=qnorm_params["num_cpus"],
                        start=qnorm_params["start"],
                        end=qnorm_params["end"])
        
        create_log_file(subcommand="quantnorm", 
                sisana_version=s_version,
                netzoopy_version=nzp_version,
                params_dict=qnorm_params, 
                filenames=outfiles)
