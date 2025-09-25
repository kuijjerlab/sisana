import yaml
import argparse
from netZooPy.panda.panda import Panda
from netZooPy.lioness.lioness import Lioness
from sisana.preprocessing import preprocess_data
from sisana.postprocessing import convert_lion_to_pickle, extract_tfs_genes
from sisana.analyze_networks import calculate_panda_degree, calculate_lioness_degree, compare_bw_groups, survival_analysis, perform_gsea, plot_volcano, plot_expression_degree, plot_heatmap, plot_clustermap, summarize
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
    A command line interface tool used to generate and analyze 
    PANDA and LIONESS networks. It works through subcommands. 
    The command 'sisana generate -p params.yaml', for example,
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

    # Add subcommands
    subparsers = parser.add_subparsers(title='Subcommands', dest='command')
    pre = subparsers.add_parser('preprocess', help='Filters expression data for parameters (e.g. genes) that are only present in at least m samples. Also filters each input file so they have the same genes and TFs across each', epilog=sisana.docs.preprocess_desc, formatter_class=argparse.RawDescriptionHelpFormatter)
    gen = subparsers.add_parser('generate', help='Generates PANDA and LIONESS networks', epilog=sisana.docs.generate_desc, formatter_class=argparse.RawDescriptionHelpFormatter)
    comb = subparsers.add_parser('combine', help='Combines indegree and outdegree files ran in batches', epilog=sisana.docs.combine_desc, formatter_class=argparse.RawDescriptionHelpFormatter)
    ext = subparsers.add_parser('extract', help='Extract edges connected to specified TFs/genes', epilog=sisana.docs.extract_desc, formatter_class=argparse.RawDescriptionHelpFormatter)
    comp = subparsers.add_parser('compare', help='Compare networks between sample groups', epilog=sisana.docs.compare_desc, formatter_class=argparse.RawDescriptionHelpFormatter)
    surv = subparsers.add_parser('survival', help='Compare survival times of individuals between sample groups', epilog=sisana.docs.survival_desc, formatter_class=argparse.RawDescriptionHelpFormatter)
    gsea = subparsers.add_parser('gsea', help='Perform gene set enrichment analysis between sample groups', epilog=sisana.docs.gsea_desc, formatter_class=argparse.RawDescriptionHelpFormatter)
    vis = subparsers.add_parser('visualize', help='Visualize the calculated degrees of each sample group', epilog=sisana.docs.visualize_desc, formatter_class=argparse.RawDescriptionHelpFormatter)
    summ = subparsers.add_parser('summarize', aliases=["summarise"], help='Summarize the outputs in an html file', epilog=sisana.docs.summarize_desc, formatter_class=argparse.RawDescriptionHelpFormatter)

    # options for preprocess subcommand
    pre.add_argument("params", type=str, help='Path to yaml file containing the parameters to use')
        
    # options for generate subcommand    
    gen.add_argument("params", type=str, help='Path to yaml file containing the parameters to use')

    # options for generate subcommand    
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
    vis.add_argument("plotchoice", type=str, choices = ["all", "quantity", "heatmap", "volcano"], nargs='?', default="all", help="The type of plot to create")   
    vis.add_argument("params", type=str, help='Path to yaml file containing the parameters to use')
    
    # options for summarize subcommand    
    summ.add_argument("logdir", nargs='?', type=str, default="./log_files/", help='Path to the directory containing the previously made log files')

    args = parser.parse_args()
      

      
    # If user wants example files, retrieve them from their installed paths
    if args.example:
        
        print("Downloading example input files from Zenodo. Please wait...")
        fetch_files()
        print("Example input files have been created in ./example_inputs/")
        sys.exit(0)
    
    if args.command != "summarize" and args.command != "summarise": 
        params = yaml.load(open(args.params), Loader=yaml.FullLoader)
    else:
        summarize(args.logdir)
        sys.exit(0)

    # Create output for temp files if one does not already exist
    os.makedirs('./tmp/', exist_ok=True)

    ########################################################
    # 1) Preprocess the data
    ########################################################

    preprocess_params = params['preprocess']
    if args.command == 'preprocess':
        
        # # Save the order of the sample names to their own file, then export the data frame without a header, since that is what is required for CLI version of PANDA
        # expdf = pd.read_csv(preprocess_params['exp_file'], sep='\t', index_col=0)
        # name_list = list(expdf.columns.values)
        
        # with open('./tmp/samples.txt', 'w') as f:
        #     for line in name_list:
        #         f.write(f"{line}\n")
        
        # Remove genes that are not expressed in at least the user-defined minimum ("number")
        results = preprocess_data(preprocess_params['exp_file'], 
                        preprocess_params['filetype'], 
                        preprocess_params['number'],
                        preprocess_params['outdir'])    
        
        fname, genes_kept, genes_removed = results[0], results[1], results[2] 
            
        removed_str = f"genes removed: {genes_removed}"
        kept_str = f"genes kept: {genes_kept}"
        
        extra_info_preprocess = [removed_str, kept_str]
        
        create_log_file("preprocess", 
            preprocess_params, 
            [fname], extra_info_preprocess)
        
    ########################################################
    # 2) Run PANDA, using the parameters from the yaml file
    ########################################################

    if args.command == 'generate':
        
        generate_params = params['generate']

        if generate_params['method'].lower() == 'panda' or generate_params['method'].lower() == 'lioness':

            # data_paths = yaml.load(open('./tmp/processed_data_paths.yml'), Loader=yaml.FullLoader)
            
            # Create output dir if one does not already exist
            panda_output_location = generate_params['pandafilepath']

            pandapath = Path(panda_output_location)
            if str(pandapath)[-4:] != ".txt":
                raise Exception("Error: Panda output file must have a .txt extension. Please edit your pandafilepath variable in your params file")
            os.makedirs(pandapath.parent, exist_ok=True)
            
            panda_obj = Panda(expression_file=generate_params['exp'], 
                motif_file=generate_params['motif'], 
                ppi_file=generate_params['ppi'], 
                save_tmp=False, 
                remove_missing=False, 
                keep_expression_matrix=True, 
                save_memory=False,
                modeProcess="intersection",
                with_header=True)

            panda_res = panda_obj.export_panda_results
            panda_res = panda_res.sort_values(by=["tf", "gene"])
            panda_res.to_csv(panda_output_location, sep=" ", index=False)
            # panda_obj.save_panda_results(panda_output_location, old_compatible=False)     
            
            print("Now calculating PANDA degrees...")
            calculate_panda_degree(inputfile=panda_output_location)
               
        if generate_params['method'].lower() == 'lioness':
            lioness_full_path = generate_params['lionessfilepath']

            if lioness_full_path[-4:] != ".npy":
                raise Exception("Error: Lioness output file must have a .npy extension. Please edit your lionessfilepath variable in your params file.")

            lionesspath_no_ext = lioness_full_path[:-4]

            if generate_params['start'] is not None:
                lionesspath_new_path = Path(f"{lionesspath_no_ext}_samples_{generate_params['start']}_to_{generate_params['end']}.npy")                
            else:
                lionesspath_new_path = Path(lioness_full_path)

            lioness_full_path = Path(lioness_full_path)

            # Run Lioness on a subset of samples if specified in the params file, otherwise run on all samples
            if generate_params['start'] is not None:
                
                
                Lioness(panda_obj, 
                           computing=generate_params['compute'], 
                           precision="double",
                           ncores=generate_params['ncores'], 
                           save_dir=lioness_full_path.parent, 
                           save_fmt="npy",
                           start=generate_params['start'],
                           end=generate_params['end'])
                        #    export_filename=f"./output/network/lioness_networks_samples_{generate_params['start']}_to_{generate_params['end']}.npy")
            else:
                Lioness(panda_obj, 
                           computing=generate_params['compute'], 
                           precision="double",
                           ncores=generate_params['ncores'], 
                           save_dir=lioness_full_path.parent, 
                           save_fmt="npy")

            # Rename the default name of the lioness output file, which is not an option of the current Lioness NetZooPy cli
            os.rename(os.path.join(lioness_full_path.parent, "lioness.npy"), lionesspath_new_path)

            #lion_loc = params['generate']['outdir'] + "lioness.npy"
            liondf = pd.DataFrame(np.load(lionesspath_new_path))            
                
            # To make the edges positive values for log2FC calculation later on, first need to transform 
            # edges by doing ln(e^w + 1), then calculate degrees. Then you can do the log2FC of degrees
            # in next step
            # 
            # This transformation is described in the paper "Regulatory Network of PD1 Signaling Is Associated 
            # with Prognosis in Glioblastoma Multiforme"
            # print("Now transforming edges...")

            # print("Datafile before transformation")
            # print(liondf.head(n=20))
            
            # lion_transformed = liondf.apply(np.vectorize(transform_edge_to_positive_val))
            
            # print("Datafile after transformation")
            # print(lion_transformed.head(n=20))        
            
            # print("LIONESS network with transformed edge values saved to " + os.path.join(params['generate']['outdir'], "lioness_transformed_edges.npy"))
            if generate_params['start'] is not None:
                pickle_path = f"./tmp/lioness_samples_{generate_params['start']}_to_{generate_params['end']}.pickle"
            else:
                pickle_path = './tmp/lioness.pickle'
                
            convert_lion_to_pickle(panda_output_location,
                                liondf,
                                "npy", 
                                './tmp/samples.txt',  
                                pickle_path,
                                start=generate_params['start'],
                                end=generate_params['end'])
            
            print("\nLIONESS networks created. Now calculating LIONESS degrees...")
            calculate_lioness_degree(inputfile=pickle_path,
                            datatype="pickle")
            print("LIONESS degrees have now been calculated.")
            
            if generate_params['start'] is not None:
                lioness_indeg_filename = f"lioness_indegree_samples_{generate_params['start']}_to_{generate_params['end']}"
                lioness_outdeg_filename = f"lioness_outdegree_samples_{generate_params['start']}_to_{generate_params['end']}"
                Path(f"./tmp/lioness_samples_{generate_params['start']}_to_{generate_params['end']}_indegree.csv").rename(f"{Path(lioness_full_path).parent}/{lioness_indeg_filename}.csv")
                Path(f"./tmp/lioness_samples_{generate_params['start']}_to_{generate_params['end']}_outdegree.csv").rename(f"{Path(lioness_full_path).parent}/{lioness_outdeg_filename}.csv")

            else:
                lioness_indeg_filename = f"lioness_indegree"
                lioness_outdeg_filename = f"lioness_outdegree"
                Path("./tmp/lioness_indegree.csv").rename(f"{Path(lioness_full_path).parent}/{lioness_indeg_filename}.csv")
                Path("./tmp/lioness_outdegree.csv").rename(f"{Path(lioness_full_path).parent}/{lioness_outdeg_filename}.csv")

            print(f"LIONESS network saved to {str(lionesspath_new_path)}")
            print(f"LIONESS degrees saved to:")
            print(f"{Path(lioness_full_path).parent}/{lioness_indeg_filename}.csv")
            print(f"{Path(lioness_full_path).parent}/{lioness_outdeg_filename}.csv")
                
        print(f"\nPANDA network saved to {panda_output_location}")
        print(f"PANDA degrees saved to:") 
        print(f"{str(panda_output_location)[:-4]}_outdegree.csv")
        print(f"{str(panda_output_location)[:-4]}_indegree.csv")
        
        outfiles = [panda_output_location,
                    f"{str(panda_output_location)[:-4]}_outdegree.csv",
                    f"{str(panda_output_location)[:-4]}_indegree.csv",
                    str(lioness_full_path),
                    f"{Path(lioness_full_path).parent}/{lioness_indeg_filename}.csv",
                    f"{Path(lioness_full_path).parent}/{lioness_outdeg_filename}.csv"]
            
        create_log_file("generate", 
                        generate_params, 
                        outfiles)
        
    ########################################################
    # 2.5) (OPTIONAL) Combine the multiple degree files into a single
    #      output file. Only used samples were "batched" when creating
    #      the single-sample networks in the previous step
    ########################################################
    if args.command == 'combine':
        
        combine_params = params['combine']
        degree_dir_path = str(Path(combine_params['degree_dir']))

        indeg_dataframes = []
        indeg_filenames = []
        outdeg_dataframes = []
        outdeg_filenames = []
        numpy_dataframes = []
        numpy_filenames = []

        with open('./tmp/samples.txt', 'r') as file:
            samplist = file.readlines()
            # Optional: Remove newline characters
            samplist = [samp.strip() for samp in samplist] 
        print(samplist)

        panda_file = pd.read_csv(combine_params['panda_file'], sep = " ")
        panda_file["edge"] = panda_file["tf"] + "-" + panda_file["gene"] 
        panda_file.index = panda_file["edge"]
    
        # with open("./tmp/samples.txt", "w") as f:
        #     for samp in expdf.columns:
        #         f.write(col + "\n")

        def _get_batched_files(df_list: list, fname_list: list, regex: str, ext: str):

            if ext == "csv":
                for file in glob.glob(f"{degree_dir_path}/{regex}"):
                    df = pd.read_csv(file, index_col=0)
                    # print(df)
                    df_list.append(df)
                    fname_list.append(file)
            else:
                for file in glob.glob(f"{degree_dir_path}/{regex}"):
                    noext = file[:-4]
                    startsamp, endsamp = int(noext.split("_")[-3]), int(noext.split("_")[-1])

                    numpy_file = np.load(file)
                    # print(len(numpy_file))
                    data = pd.DataFrame(numpy_file)

                    data.columns = samplist[startsamp-1:endsamp]
                    data.index = panda_file.index
          
                    # print(panda_file)
                    # print(data)
                    df_list.append(data)
                    fname_list.append(file)

        # Combine the degree files automatically, since they are relatively small.
        # Combining networks may run into memory issues, so it's optional
        _get_batched_files(outdeg_dataframes, outdeg_filenames, "lioness_outdegree_samples_*_to_*.csv", "csv")
        _get_batched_files(indeg_dataframes, indeg_filenames, "lioness_indegree_samples_*_to_*.csv", "csv")    
        
        combined_indeg = pd.concat(indeg_dataframes, axis=1)
        combined_outdeg = pd.concat(outdeg_dataframes, axis=1)
        
        combined_indeg.to_csv(f"{degree_dir_path}/lioness_indegree.csv", index=True)
        combined_outdeg.to_csv(f"{degree_dir_path}/lioness_outdegree.csv", index=True)

        print(f"File created: {degree_dir_path}/lioness_indegree.csv")
        print(f"File created: {degree_dir_path}/lioness_outdegree.csv")
                    
        if combine_params["delete_intermediate_files"] == True:
            [os.remove(file) for file in indeg_filenames]
            [os.remove(file) for file in outdeg_filenames]
                    
        if params["combine"]["networks"] == True:
            _get_batched_files(numpy_dataframes, numpy_filenames, "lioness_networks_samples_*_to_*.npy", "npy")                  
            combined_nw = pd.concat(numpy_dataframes, axis=1)
            
            # print(combined_nw)
            
            pickle_path = './tmp/lioness.pickle'
            np_path = f"{degree_dir_path}/lioness_network.npy"
            
            with open("./tmp/combined_samples.txt", "w") as f:
                for col in combined_nw.columns:
                    f.write(col + "\n")

            convert_lion_to_pickle(combine_params['panda_file'],
                        combined_nw,
                        "npy", 
                        './tmp/combined_samples.txt',  
                        pickle_path)
            
            combined_nw = combined_nw.to_numpy()
            np.save(np_path, combined_nw)
            
            # combined_nw.to_csv(np_path, index=True)
            print(f"File created: {np_path}")
                
            if combine_params["delete_intermediate_files"] == True:
                [os.remove(file) for file in numpy_filenames]
                        

    ########################################################
    # 3) Compare degree (or expression) between sample groups
    ########################################################
        
    if args.command == "compare":     
        compare_means_params = params['compare']

        outfiles = compare_bw_groups(datafile=compare_means_params["datafile"], 
                                    mapfile=compare_means_params["mapfile"], 
                                    datatype=compare_means_params["datatype"], 
                                    groups=compare_means_params["groups"],
                                    testtype=compare_means_params["testtype"], 
                                    filetype=compare_means_params["filetype"],
                                    rankby_col=compare_means_params["rankby"],
                                    outdir=compare_means_params["outdir"])
        
        create_log_file("compare_means", 
            compare_means_params, 
            outfiles)
    
    ########################################################
    # 4) Perform gene set enrichment analysis
    ########################################################   
        
    if args.command == 'gsea':    
        gsea_params = params["gsea"]
        
        outfiles = perform_gsea(genefile=gsea_params["genefile"], 
                        gmtfile=gsea_params["gmtfile"], 
                        geneset=gsea_params["geneset"], 
                        outdir=gsea_params["outdir"])
        
        create_log_file("gsea", 
            gsea_params, 
            outfiles)
    
    ########################################################
    # 5) Visualize results
    ########################################################       

    if args.command == "visualize":                  

        if args.plotchoice == "volcano":    
            volcano_params = params["visualize"]["volcano"]

            outfiles, down_gene_count, up_gene_count = plot_volcano(statsfile=volcano_params["statsfile"],
                         diffcol=volcano_params["diffcol"],
                         adjpcol=volcano_params["adjpcol"],
                         adjpvalthreshold=volcano_params["adjpvalthreshold"],
                         xaxisthreshold=volcano_params["xaxisthreshold"],
                         difftype=volcano_params["difftype"],
                         genelist=volcano_params["genelist"],
                         outdir=volcano_params["outdir"],
                         top=False)      
            
            down_gene_str = f"number of genes up in group 1: {down_gene_count}"
            up_gene_str = f"number of genes up in group 2: {up_gene_count}"
            
            extra_info_num_genes = [down_gene_str, up_gene_str]
            
            create_log_file("volcano_plot", 
                volcano_params, 
                [outfiles], extra_info_num_genes)
    
        if args.plotchoice == "quantity":   
            quantity_params = params["visualize"]["quantity"]
            
            if quantity_params["genelist"] != None:
                outfiles = plot_expression_degree(datafile=quantity_params["datafile"],
                            filetype=quantity_params["filetype"], 
                            statsfile=quantity_params["statsfile"], 
                            metadata=quantity_params["metadata"],
                            plottype=quantity_params["plottype"],
                            groups=quantity_params["groups"],
                            colors=quantity_params["colors"],
                            prefix=quantity_params["prefix"],
                            yaxisname=quantity_params["yaxisname"],
                            outdir=quantity_params["outdir"],
                            genelist=quantity_params["genelist"],
                            top=False)   
            else:
                outfiles = plot_expression_degree(datafile=quantity_params["datafile"],
                            filetype=quantity_params["filetype"], 
                            statsfile=quantity_params["statsfile"], 
                            metadata=quantity_params["metadata"],
                            plottype=quantity_params["plottype"],
                            groups=quantity_params["groups"],
                            colors=quantity_params["colors"],
                            prefix=quantity_params["prefix"],
                            yaxisname=quantity_params["yaxisname"],
                            outdir=quantity_params["outdir"],
                            genelist=quantity_params["genelist"],
                            numgenes=quantity_params["numgenes"],
                            top=True)   
                
            create_log_file("quantity_plot", 
                quantity_params, 
                [outfiles])               
                
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
            heatmap_params = params["visualize"]["heatmap"]

            outfiles = plot_clustermap(datafile=heatmap_params["datafile"],
                        filetype=heatmap_params["filetype"], 
                        metadata=heatmap_params["metadata"],
                        genelist=heatmap_params["genelist"],
                        column_cluster=heatmap_params["column_cluster"],
                        row_cluster=heatmap_params["row_cluster"],
                        prefix=heatmap_params["prefix"],
                        outdir=heatmap_params["outdir"],
                        plot_gene_names=heatmap_params["plot_gene_names"],
                        plot_sample_names=heatmap_params["plot_sample_names"],
                        category_label_columns=heatmap_params["category_label_columns"],
                        category_column_colors=heatmap_params["category_column_colors"],                       
                        top=False)   
            
            create_log_file("heatmap", 
                heatmap_params, 
                outfiles)  
            
    ########################################################
    # (Optional) Extract edges that connect to specific TFs/genes
    ########################################################

    if args.command == 'extract':
        extract_params = params["extract"]

        outfiles = extract_tfs_genes(pickle=extract_params["pickle"], 
                         datatype=args.extractchoice, 
                         sampnames=extract_params["sampnames"],
                         symbols=extract_params["symbols"], 
                         outdir=extract_params["outdir"])
        
        create_log_file("extract", 
                extract_params, 
                [outfiles])  
        
    ########################################################
    # (Optional) Perform survival analysis
    ########################################################
   
    if args.command == "survival":     
        compare_survival_params = params['survival']

        try:
            outfiles = survival_analysis(metadata=compare_survival_params["metadata"],
                            filetype=compare_survival_params["filetype"], 
                            sampgroup_colname=compare_survival_params["sampgroup_colname"],
                            alivestatus_colname=compare_survival_params["alivestatus_colname"],
                            days_colname=compare_survival_params["days_colname"],
                            groups=compare_survival_params["groups"],
                            outdir=compare_survival_params["outdir"],
                            appendname=compare_survival_params["appendname"])
        except:
            outfiles = survival_analysis(metadata=compare_survival_params["metadata"],
                            filetype=compare_survival_params["filetype"], 
                            sampgroup_colname=compare_survival_params["sampgroup_colname"],
                            alivestatus_colname=compare_survival_params["alivestatus_colname"],
                            days_colname=compare_survival_params["days_colname"],
                            groups=compare_survival_params["groups"],
                            outdir=compare_survival_params["outdir"])
        fnames, pval, sig = outfiles[0], outfiles[1], outfiles[2] 
        
        pval_str = f"p-value: {pval}"
        sig_str = f"significant?: {sig}"
        
        extra_info = []
        extra_info.append(pval_str)
        extra_info.append(sig_str)
        
        create_log_file("compare_survival", 
            compare_survival_params, 
            [fnames], extra_info)
