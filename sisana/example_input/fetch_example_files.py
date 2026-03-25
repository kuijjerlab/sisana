import requests
import os
from tqdm import tqdm
import sys

def fetch_files():
    """
    Description:
        This code fetches the example input files from Zenodo and saves them in a directory called ./example_inputs/.
        
    Parameters:
    -----------
        - None
    
    Returns:
    -----------
        - Nothing
    """
    
    os.makedirs('./example_inputs/', exist_ok=True)
    
    # zenodo tag (must update with each new release). The tag is the XXXXXXXX in the following: https://zenodo.org/records/XXXXXXXX/files/*
    tag = "17190642"
    
    sisana_suffixes = ['/files/BRCA_TCGA_200_LumA_LumB_samps_mapping_w_header.csv',
            '/files/BRCA_TCGA_200_LumA_LumB_samps_survival_data.csv',
            '/files/BRCA_TCGA_20_LumA_LumB_samps_5000_genes_exp.tsv',
            '/files/c2.cp.kegg_medicus.v2023.2.Hs.symbols.gmt',
            '/files/c2.cp.reactome.v2023.2.Hs.symbols.gmt',
            '/files/clustering_template.ipynb',
            '/files/genes_to_extract.txt',
            '/files/Hallmark.v2023.2.Hs.symbols.gmt',
            '/files/heatmap_genes.txt',
            '/files/lioness_df_indegree_3_decimal_places_subset_200_LumALumB_samps.csv',
            '/files/params.yml',
            '/files/params.yml',
            '/files/quantity_plot_genes.txt',
            '/files/volcano_plot_genes.txt']
            
    ### From the SPONGE Zenodo repo:
    sponge_suffixes = ['/files/ppi_prior_2024.tsv',
                       '/files/motif_prior_names_2024.tsv']

    sisana_record_link = 'https://zenodo.org/records/17190642'
    sponge_record_link = 'https://zenodo.org/records/13628785'
    # file_suffix = '/files/BRCA_TCGA_200_LumA_LumB_samps_mapping_w_header.csv'

    # Redirecting only works when accessing the main site
    sisana_zenodo_url = requests.get(sisana_record_link, allow_redirects=True)
    sponge_zenodo_url = requests.get(sponge_record_link, allow_redirects=True)
    # Now use the redirected url to access files
    # redirected_zenodo_url = requests.get(zenodo_url.url + file_suffix)
    # print(redirected_zenodo_url.text)
    # sys.exit(0)

    def _download_files(suffixes: list, zenodo_url: str, method_name: str):
        curfile = 1
        for suffix in suffixes:
            filename = suffix.split("/")[-1]
            print(f"Downloading file {curfile} of {len(suffixes)} ({filename}) from the {method_name} Zenodo repository...")
            
            print(zenodo_url.url + suffix)
            redirected_zenodo_url = requests.get(zenodo_url.url + suffix, stream=True, allow_redirects=True)
            # r = requests.get(redirected_zenodo_url, stream=True, allow_redirects=True)
            # fname = os.path.basename(redirected_zenodo_url)
            
            # Sizes in bytes.
            total_size = int(redirected_zenodo_url.headers.get("content-length", 0))
            block_size = 1024       
            
            with tqdm(total=total_size, unit="B", unit_scale=True) as progress_bar:
                with open(f'./example_inputs/{filename}', 'wb') as f:
                    for chunk in redirected_zenodo_url.iter_content(block_size):
                        progress_bar.update(len(chunk))
                        f.write(chunk)
            print("\n")
            curfile += 1
            f.close()

    _download_files(sisana_suffixes, sisana_zenodo_url, "SiSaNA")
    _download_files(sponge_suffixes, sponge_zenodo_url, "SPONGE")

    # for suffix in sisana_suffixes:
    #     filename = suffix.split("/")[-1]
    #     print(f"Downloading file {curfile} of {len(sisana_suffixes)}: {filename}")
        
    #     print(zenodo_url.url + suffix)
    #     redirected_zenodo_url = requests.get(zenodo_url.url + suffix)
    #     print(redirected_zenodo_url.text)
    #     # r = requests.get(redirected_zenodo_url, stream=True, allow_redirects=True)
    #     # fname = os.path.basename(redirected_zenodo_url)
        
    #     # Sizes in bytes.
    #     total_size = int(redirected_zenodo_url.headers.get("content-length", 0))
    #     block_size = 1024       
        
    #     with tqdm(total=total_size, unit="B", unit_scale=True) as progress_bar:
    #         with open(f'./example_inputs/{filename}', 'wb') as f:
    #             for chunk in redirected_zenodo_url.iter_content(block_size):
    #                 progress_bar.update(len(chunk))
    #                 f.write(chunk)
    #     print("\n")
    #     curfile += 1
    #     f.close()