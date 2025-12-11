# SiSaNA - Single Sample Network Analysis

<p align="center">  
  <img src="https://github.com/kuijjerlab/sisana/blob/main/docs/sisana_logo.png" width="300" />
</p>

## Read about SiSaNA
[SiSaNA can currently be found on bioRxiv](https://www.biorxiv.org/content/10.1101/2025.11.06.680212v1)

## Table of contents
- [SiSaNA introduction](#sisana-introduction)
  - [What is SiSaNA?](#what-is-sisana)
  - [What are PANDA and LIONESS?](#what-are-panda-and-lioness)
  - [How do you interpret the gene regulatory networks?](#how-do-you-interpret-the-gene-regulatory-networks)
- [Requirements](#requirements-for-using-sisana)
- [Installation](#installation)
- [Before you begin](#before-you-begin)
  - [Example input files](#example-input-files)
  - [Help documentation](#sisana-help-documentation)
  - [Setting up your params.yml file](#setting-up-your-paramsyml-file)
- [Running SiSaNA](#running-sisana)
  - [Step 1: Pre-processing of data](#step-1-pre-processing-of-data)
  - [Step 2: Reconstruction and analysis of networks](#step-2-reconstruction-and-analysis-of-networks)
  - [Step 3: Comparing two experimental groups](#step-3-comparing-two-experimental-groups)
  - [Step 4: Survival analysis](#step-4-survival-analysis)
  - [Step 5: Performing gene set enrichment analysis](#step-5-performing-gene-set-enrichment-analysis)
  - [Step 6: Visualization of results](#step-6-visualization-of-results)
  - [Step 7: Summarize your results](#step-7-summarize-your-results)

# SiSaNA introduction

## What is SiSaNA? 
SiSaNA is a command line tool that utiliizes the PANDA and LIONESS algorithms from the netZooPy module to reconstruct single sample regulatory networks. Using SiSaNA, users can easily calculate in- and out-degree for each of the reconstructed networks. Additionally, SiSaNA can compare the expression/degree between groups of interest, including performing statistical tests, visualizing the results (volcano plots, boxplots, violin plots, and heatmaps), and compare the survival between groups. All this is accomplished via the command line, with little to no prior programming experience required. 

## What are PANDA and LIONESS?
PANDA is a tool developed to reconstruct gene regulatory networks from bulk input data (such as RNA-Seq data). It uses a message passing approach, along with prior protein-protein interaction data and prior regulatory knowledge information to refine a single network that reflects the regulatory landscape of the input samples. LIONESS utilizes the PANDA algorithm and iteratively removes one sample (with replacement), then uses PANDA to reconstruct a network with all samples minus one. It then uses differences in the original PANDA network and the newly created network with n-1 samples to reconstruct single-sample regulatory networks.

<p align="center">  
  <img src="https://github.com/kuijjerlab/sisana/blob/main/docs/panda_lioness_example.png" width="500" />
</p>

## How do you interpret the gene regulatory networks?
These networks consist of two types of nodes: transcription factors (TFs) and genes, with an edge that connects each TF to each gene. The weight (or value) of the edge denotes the likelihood of a TF to regulate that gene. A larger edge weight means a higher likelihood of regulation. 

SiSaNA will also calculate in-degree and out-degree, where in-degree is the sum of edge weights coming in to a gene while out-degree is the sum of edge weights coming out from a TF. As you may have guessed, a larger in-degree means that, on average, a gene is being more highly regulated in the modeled disease. Meanwhile, a larger out-degree means that, on average, that TF is regulating more genes in the modeled disease.

<p align="center">  
  <img src="https://github.com/kuijjerlab/sisana/blob/main/docs/network_example.png" width="300" />
</p>

## Requirements for using SiSaNA
 - python v3.9.19 (see installation steps for creating a conda environment with this specific Python version). SiSaNA should work with versions of Python 3.9.0 or greater, but as it has been written and tested on this version, we will use 3.9.19.
   
## Installation

1. Create a conda virtual environment with python version 3.9.19. Note: You need to substitute the path you want on your own system for the --prefix argument
```
conda create --prefix </path/to/env-name> python=3.9.19
```

2. Enter the conda environment
```
conda activate </path/to/env-name>
```

3. Install SiSaNA via the pip package installer
```
pip3 install sisana
```

4. Create a directory for the analysis and move into the analysis directory
```
mkdir sisana
cd sisana
```

# Before you begin

## Pipeline overview 
<p align="center">  
  <img src="https://github.com/kuijjerlab/sisana/blob/main/docs/sisana_pipeline_overview_v12.png" width="400" />
</p>

## Example input files
Example input files can be obtained using the command
```
sisana -e
```
These files will be downloaded from [Zenodo](https://zenodo.org/records/17190643) and stored in a directory called "example_inputs". One of these example files is the params.yml file, which can be used as a template and edited for your own data (see next section). Each user-defined parameter in the params.yml file is documented with a comment to explain the function of the parameter. The comments do not need to be removed prior to running SiSaNA. The files in this example_inputs directory can be used in the commands listed down below.

#### ❓Question 1: What is the structure of each of the three input files for SiSaNA's preprocessing step? What do the columns represent in each of the files? 

#### ❓Question 2: How many genes are being used as input in the expression file? How many unique TFs are in the prior motif file? How many unique genes in the prior motif file?

## SiSaNA help documentation
To view help documentation on which subcommands are available, the following can be used:
```
sisana -h
```

For further information on these subcommands, simply put the name of the subcommand before the `-h`
```
sisana <subcommand> -h
```

## Setting up your params.yml file
The most important thing to get right in order to correctly run SiSaNA is the structure of your params.yml file. SiSaNA comes with a params.yml file that is annotated to explain the function of each argument. The params.yml file is separated into 'chunks' that reflect the same subcommands available in SiSaNA on the command line. For each step of SiSaNA, you will need to use the correct subcommand, as well as have the parameters set up in the params.yml file.

In the below example, the user is running the "preprocess" step of SiSaNA. They have specified the paths to the input files as well as the value for the number of samples a gene must be expressed in (in their case, 5), along with the path to the output directory in which to store their results.
![Pipeline overview](docs/params_example_v2.png)

# Running SiSaNA

## Step 1: Pre-processing of data
The "preprocess" subcommand is the first stage of SiSaNA, where it preprocess the input data to get it in a format that the PANDA and LIONESS algorithms can handle. This will likely involve the removal of genes or transcription factors that are not consistent across files. Information regarding the removal of these factors is given at the end of the preprocessing step.

#### Example command
```
sisana preprocess ./example_inputs/params.yml
```

#### Outputs
Three files, one for each of the three filtered input files. 
<br />
<br />


## Step 2: Reconstruction and analysis of networks
This second SiSaNA stage, "generate", uses the PANDA and LIONESS algorithms of netZooPy to reconstruct gene regulatory networks. Documentation for netZooPy can be found at https://github.com/netZoo/netZooPy/tree/master. It then performs basic analyses of these networks by calculating in-degree of genes (also called gene targeting scores) and out-degree of transcription factors (TFs).

#### Example command
```
sisana generate ./example_inputs/params.yml
```

#### Outputs
1. lioness.npy, which contains all calculated edges for each sample
2. lioness.pickle, which is the same thing, just serialized to make reading into python quicker
3. A file containing the calculated indegree and another file with the outdegree of each gene and transcription factor, respectively.
<br />
<br />

#### ❓Question 3: After running the generate step, a few files are automatically generated. Which of them might you use to determine which genes were being most highly regulated?

## Step 3: Comparing two experimental groups
The next stage in SiSaNA, "compare", is used to find out how the in-degree and out-degree differ between each group. SiSaNA offers multiple ways to do this comparison, including t-tests (and Mann-Whitney tests) and paired t-tests (and Wilcoxon paired t-tests).

#### Example commands
To compare the values between two groups in order to identify differentially expressed genes or differential degrees, you can use the following command:
```
sisana compare ./example_inputs/params.yml
```
<br />
<br />

#### ❓(Ignore this question...) ~~Question 4: Following the comparison of the expression of the two Luminal breast cancer groups (the default setup in the params.yml file), take a look at the output .txt file. What are some of the most differentially expressed genes?~~

#### ❓Question 5: What are some of the top genes that differ in the amount they are being regulated? 

## Step 4: Survival analysis
For performing survival analyses, you can use a command like this:
```
sisana survival ./example_inputs/params.yml
```
<br />
<p align="center">
  <img src="https://github.com/kuijjerlab/sisana/blob/main/docs/LumA_v_LumB_survival_plot.png" width="500" />
</p>
<br />
<br />

#### ❓Question 6: . Take a look at the output of the "sisana survival" command. Do you know how to interpret this image? Additionally, there is a p-value listed there. Do you think we can conclude based on this p-value that there is a significant difference in the survival of Luminal A vs Luminal B patients?



## Step 5: Performing gene set enrichment analysis 
"sisana gsea" is used to perform gene set enrichment analysis (GSEA) to identify pathways that are differentially regulated based on the gene targeting scores. It uses the ranks of genes found in the previous step (sisana compare) as input.

#### Example commands
```
sisana gsea ./example_inputs/params.yml
```
<br />

<p align="center">
  <img src="https://github.com/kuijjerlab/sisana/blob/main/docs/comparison_mw_between_LumA_LumB_degree_ranked_mediandiff_GSEA_Hallmarks_basic_enrichment_dotplot.png" width="500" />
</p>
<br />
<br />

#### ❓Question 7: Following running the "sisana gsea" step, what are the three top pathways? Run this step again, but with the Reactome gene set instead of the Hallmarks one.

## Step 6: Visualization of results
The "visualize" command allows you to visualize the results of your analysis on publication-ready figures. There are multiple types of visualization you can perform, including generating volcano plots...
```
sisana visualize volcano ./example_inputs/params.yml
```
<p align="center">
  <img src="https://github.com/kuijjerlab/sisana/blob/main/docs/volcano_plot_adjp_0.25.png" width="500" />
</p>

#### ❓Question 8. Now let's look at the output of the "sisana visualize volcano" command. What do the dots on the left and right side of the plot represent? Do you recognize any of these genes?

...making boxplots or violin plots of expression/degrees...
```
sisana visualize quantity ./example_inputs/params.yml
```
<<<<<<< HEAD

#### ❓Question 9. Take 5 of the top differential indegrees (which you calculated in step 4) and make a box plot for those genes.
=======
<p align="center">
  <img src="https://github.com/kuijjerlab/sisana/blob/main/docs/LumA_LumB_indegree_box_plot.png" width="500" />
</p>
>>>>>>> parent of 853eb82 (Removed images for OBiWoW workshop)

...and creating heatmaps
```
sisana visualize heatmap ./example_inputs/params.yml
```
<p align="center">
  <img src="https://github.com/kuijjerlab/sisana/blob/main/docs/TCGA_200_LumA_LumB_samps_heatmap.png" width="500" />
</p> 
<br />
<br />

## Step 7: Summarize your results
The final stage of SiSaNA, "summarize", takes all the created images and outputs them in a single html file for convenience. This can then be opened in a web browser. Please note that you must be in the directory containing the log_files subdirectory for this command to work.
```
sisana summarize
```
<p align="center">
  <img src="https://github.com/kuijjerlab/sisana/blob/main/docs/example_html_output.png" width="750" />
</p> 

#### ❓Question 10. How many TFs have edges that connect to (regulate) the GDI2 gene in your network? (Hint: Use the help function of SiSaNA (sisana -h) to find a subcommand that can help with this).
