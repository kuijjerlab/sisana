
def get_default_params():
    """
    Description:
        This function generates a dictionary that holds the default parameters for SiSaNA. These values then get overwritten
        by the values the user has in the params.yml file as needed.
        
    Parameters:
    -----------
        - None
        
    Returns:
    -----------
        - dict: dictionary in the form of the params.yml input, where keys are stages of SiSaNA and values are the parameters 
                set for each stage (except for visualize, which is a nested dict at this time)
    """
    def_params = {}
    
    def_params["preprocess"] = {}
    def_preprocess = def_params["preprocess"]
    def_preprocess["number"] = 5
    def_preprocess["outdir"] = "./output/preprocess"
    
    def_params["generate"] = {}
    def_generate = def_params["generate"]
    def_generate["method"] = "lioness"
    def_generate["modeProcess"] = "intersection"
    def_generate["pandafilepath"] = "./output/network/panda_network.txt"
    def_generate["compute"] = "cpu"
    def_generate["ncores"] = 2
    def_generate["lionessfilepath"] = "./output/network/lioness_networks.npy"
    def_generate["start"] = None
    def_generate["end"] = None
    
    def_params["compare"] = {}
    def_compare = def_params["compare"]
    def_compare["datatype"] = "degree"
    def_compare["testtype"] = "mw"
    def_compare["rankby"] = "mediandiff"
    def_compare["outdir"] = "./output/compare_means/"
    
    def_params["survival"] = {}
    def_survival = def_params["survival"]  
    def_survival["outdir"] = "./output/survival/"  

    def_params["gsea"] = {}
    def_survival = def_params["gsea"]  
    def_survival["outdir"] = "./output/gsea/"  

    def_params["visualize"] = {}
    def_params["visualize"]["volcano"] = {}
    def_volcano = def_params["visualize"]["volcano"]
    def_volcano["adjpcol"] = "FDR"
    def_volcano["xaxisthreshold"] = 50
    def_volcano["adjpvalthreshold"] = 0.05
    def_volcano["difftype"] = "median"
    def_volcano["outdir"] = "./output/volcano/"
    def_volcano["top"] = True 
    def_volcano["numlabels"] = 15

    def_params["visualize"]["quantity"] = {}
    def_quantity = def_params["visualize"]["quantity"]
    def_quantity["plottype"] = "boxplot"
    def_quantity["colors"] = ["cornflowerblue", "orange"]
    def_quantity["outdir"] = "./output/plot_quantity/"
    def_quantity["prefix"] = "results"
    def_quantity["genelist"] = None
    def_quantity["numgenes"] = 5

    def_params["visualize"]["heatmap"] = {}
    def_heatmap = def_params["visualize"]["heatmap"]
    def_heatmap["column_cluster"] = False
    def_heatmap["row_cluster"] = True
    def_heatmap["plot_gene_names"] = True
    def_heatmap["plot_sample_names"] = False
    def_heatmap["outdir"] = "./output/heatmap/"
    def_heatmap["prefix"] = "results"
    # def_heatmap["top"] = True    

    def_params["extract"] = {}
    def_extract = def_params["extract"]  
    def_extract["pickle"] = "./tmp/lioness.pickle"  
    def_extract["sampnames"] = "./tmp/samples.txt"  
    def_extract["outdir"] = "./output/extract/"  

    def_params["combine"] = {}
    def_combine = def_params["combine"]  
    def_combine["degree_dir"] = "./output/network/"  
    def_combine["panda_file"] = "./output/network/panda_network.txt"  
    def_combine["networks"] = True
    def_combine["delete_intermediate_files"] = True

    return def_params