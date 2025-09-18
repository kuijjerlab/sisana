import base64
import os
import sys
import yaml
from pathlib import Path

def summarize(logdir: str):
    """
    Description:
        Uses the log files generated to find and summarize the resulting figures in html format
        
    Parameters:
    -----------
        - logdir: str, Path to the directory that contains the log files 

    Returns:
    -----------
        - None
    """
    
    def _read_log_files(logdir: str):
        """
        Description:
            Parses through the files in the log_files directory to identify created image files
            
        Parameters:
        -----------
            - logdir: str, Path to the directory that contains the log files 

        Returns:
        -----------
            - dictionary of image files found, keyed on their analysis type
        """
        
        def __populate_analysis_info_dict(yaml_file: dict, analysis_info_dict: dict, analysis_type: str):
            """
            Description:
                Creates a dictionary keyed on analysis type (preprocess, generate, etc.) and with
                values of the input, output, and additional info (if additional info was returned)
                by the analysis function previously. Note that the structure of this created dict
                is similar to the input yaml_file, but the input does not contain the keys based on
                the analysis type
                
            Parameters:
            -----------
                - yaml_file: str, yaml object that includes the read-in log yaml file
                - analysis_info_dict: dict, keyed on analysis type (preprocess, generate, etc.) and with
                  values of the input, output, and additional info
                - analysis_type: str, 

            Returns:
            -----------
                - dictionary keyed on analysis type (preprocess, generate, etc.) and with values of the 
                  input, output, and additional info
            """
            analysis_info_dict[analysis_type] = {}
            analysis_info_dict[analysis_type]["input"] = yaml_file["Parameters used"]

            # Check if the object should be processed as a string or list, then fill in the 
            # "input" and "output" keys of that analysis type
            if not isinstance(log_yaml["Files generated"], str):
                analysis_info_dict[analysis_type]["output"] = [outfile for outfile in log_yaml["Files generated"] if outfile.endswith("png")]
            else:   
                if log_yaml["Files generated"].endswith("png"):
                    analysis_info_dict[analysis_type]["output"] = yaml_file["Files generated"]
        
            # "Additional information" includes additional relevant information that would be 
            # useful for reference for the user, including some basic statistics that are not
            # normally part of the main output of the analysis, such as the number of genes that
            # are retained after filtering for genes not expressed across X samples
            try:
                analysis_info_dict[analysis_type]["Additional information"] = yaml_file["Additional information"]
            except:
                analysis_info_dict[analysis_type]["Additional information"] = "None"
        
            try:
                if not isinstance(log_yaml["Additional information"], str):
                    for k,v in log_yaml["Additional information"].items():
                        analysis_info_dict[analysis_type]["Additional information"][k] = v
                else:   
                    analysis_info_dict[analysis_type]["Additional information"] = log_yaml["Additional information"]
            except Exception:
                pass
        
            return(analysis_info_dict)

        log_files = os.listdir(logdir)
        analysis_info_dict = {}

        for i in log_files:
            
            # Read in yaml file for each log file and store into the analysis_info_dict
            filepath = os.path.join(logdir, i)
            log_yaml = yaml.load(open(filepath), Loader=yaml.FullLoader)
            
            if i == "compare_survival_log.txt":
                __populate_analysis_info_dict(log_yaml, analysis_info_dict, "survival")
                
            if i == "quantity_plot_log.txt":
                __populate_analysis_info_dict(log_yaml, analysis_info_dict, "quantity")
                
            elif i == "heatmap_log.txt":
                __populate_analysis_info_dict(log_yaml, analysis_info_dict, "heatmap")
                
            elif i == "volcano_plot_log.txt":
                __populate_analysis_info_dict(log_yaml, analysis_info_dict, "volcano")
            
            elif i == "gsea_log.txt":
                __populate_analysis_info_dict(log_yaml, analysis_info_dict, "gsea")

            elif i == "preprocess_log.txt":
                __populate_analysis_info_dict(log_yaml, analysis_info_dict, "preprocess")

            else:
                pass

        return(analysis_info_dict)

    def _create_img_tag(analysis_type_dict):
        """
        Description:
            Creates the image tag in html format for each analysis type
        """
        images = []
        
        if not isinstance(analysis_type_dict["output"], str):
            for i in analysis_type_dict["output"]:
                data_uri = base64.b64encode(open(i, 'rb').read()).decode('utf-8')
                img_tag = '<center><img src="data:image/png;base64,{0}" style="width:600px;height:auto;" /></center>'.format(data_uri)
                images.append(img_tag)
        else:
            data_uri = base64.b64encode(open(analysis_type_dict["output"], 'rb').read()).decode('utf-8')
            img_tag = '<center><img src="data:image/png;base64,{0}" style="width:600px;height:auto;" /></center>'.format(data_uri)
            images.append(img_tag)
        return('\n'.join(images))
    
    ################################
    # Main driver code
    ################################
    all_analyses_info_dict = _read_log_files(logdir)
    
    # Create the html contents for each analysis type. If that analysis log yaml file does not exist, then the analysis type will
    # not be in the all_analyses_info_dict, so it is left blank
    if 'preprocess' in all_analyses_info_dict: 
        preprocess_add_info = {key: val for i in all_analyses_info_dict["preprocess"]["Additional information"] for key, val in i.items()}    

        preprocess_results = f"""
            <div>
                <h2>Preprocessing</h2>
                <p>After removing <strong>{preprocess_add_info["genes removed"]}</strong> genes that were not present in at least {all_analyses_info_dict["preprocess"]["input"]["number"]} samples, <strong>{preprocess_add_info["genes kept"]}</strong> genes were retained for the network reconstruction step.</p>
            </div>
        """  
    else: preprocess_results = ""
    
    if 'survival' in all_analyses_info_dict: 
        survival_add_info = {key: val for i in all_analyses_info_dict["survival"]["Additional information"] for key, val in i.items()}    
        survival_results = f"""
        <div>
            <h2>Survival plot</h2>
            <p>Your groups differ in survival with a p-value of <strong>{survival_add_info["p-value"]:.3f}</strong>. This indicates that there is a <strong>{survival_add_info["significant?"]}</strong> difference in the survival between the groups.</p>
            {_create_img_tag(all_analyses_info_dict["survival"])}
        </div>
        """   
    else: survival_results = ""

    if 'heatmap' in all_analyses_info_dict: 
        
        if all_analyses_info_dict["heatmap"]["input"]["column_cluster"] == False:
            colclust_str = "Column clustering was not performed."
        else:
            colclust_str = "Column clustering was performed."
            
        if all_analyses_info_dict["heatmap"]["input"]["row_cluster"] == False:
            rowclust_str = "Row clustering was not performed."
        else:
            rowclust_str = "Row clustering was performed."            
            
        heatmap_results = f"""
            <div>
                <h2>Quantity plot</h2>
                <p>Below you will find the heatmap created for visualizing the genes. {colclust_str} {rowclust_str}</p>
                {_create_img_tag(all_analyses_info_dict["heatmap"])}
            </div>
        """  
    else: heatmap_results = ""    
    

    if 'quantity' in all_analyses_info_dict: 
        if all_analyses_info_dict["quantity"]["input"]["numgenes"] is not None:
            ngenes_string = f"The {all_analyses_info_dict['quantity']['input']['numgenes']} top genes were plotted."
        else:
            ngenes_string = ""
            
        quantity_results = f"""
            <div>
                <h2>Quantity plot</h2>
                <p>Below you will find the {all_analyses_info_dict["quantity"]["input"]["plottype"]} created for visualizing the {all_analyses_info_dict["quantity"]["input"]["yaxisname"]}. {ngenes_string}</p>
                {_create_img_tag(all_analyses_info_dict["quantity"])}
            </div>
        """  
    else: quantity_results = ""
        
    if 'volcano' in all_analyses_info_dict: 
        volcano_add_info = {key: val for i in all_analyses_info_dict["volcano"]["Additional information"] for key, val in i.items()}    

        volcano_results = f"""
            <div>
                <h2>Volcano plot</h2>
                <p>To interpret this plot, pay attention to the TFs/genes that are colored. These are the genes that are below the <strong>{all_analyses_info_dict["volcano"]["input"]["adjpvalthreshold"]}</strong> FDR threshold and greater (absolute value) than the threshold of <strong>{all_analyses_info_dict["volcano"]["input"]["xaxisthreshold"]}</strong> set for the difference in {all_analyses_info_dict["volcano"]["input"]["difftype"]} degree. These may be genes that are important in distinguishing your two groups from one another.</p>
                <p>There are <strong>{volcano_add_info["number of genes up in group 1"]}</strong> higher in LumA and <strong>{volcano_add_info["number of genes up in group 2"]}</strong> genes higher in LumB.</p>

                {_create_img_tag(all_analyses_info_dict["volcano"])}
            </div>
        """ 
    else: volcano_results = ""
        
    if 'gsea' in all_analyses_info_dict: 
        gsea_results = f"""
        <div>
            <h2>GSEA results</h2>
            <p>Below you will find the top pathways that are enriched between the two sample groups. The GMT file used for generating these results was {Path(all_analyses_info_dict["gsea"]["input"]["gmtfile"]).name}.</p>
            {_create_img_tag(all_analyses_info_dict["gsea"])}
        </div>
        """   
    else: gsea_results = ""
       

    html_content = f"""
        <!DOCTYPE html>
        <html lang="en">
        <head>
            <meta charset="UTF-8">
            <meta name="viewport" content="width=device-width, initial-scale=1.0">
            <title>Python Results</title>
            <style>
                body {{ font-family: sans-serif; margin: 20px; }}
                h1 {{ color: #333; }}
                .result-item {{ background-color: #f0f0f0; padding: 10px; margin-bottom: 5px; border-radius: 5px; }}
            </style>
        </head>
        <body>
            <h1>SiSaNA results</h1>
            {preprocess_results}
            {survival_results}
            {quantity_results}
            {volcano_results}
            {heatmap_results}
            {gsea_results}
        </body>
        </html>
        """
        
    file_name = "results_summarized.html"
    with open(file_name, "w", encoding="utf-8") as f:
        f.write(html_content)
        
    print(f"Results have been summarized in {file_name}")