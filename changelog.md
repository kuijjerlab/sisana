### 1.5.3
- Fixed a bug where the visualization commands would fail due to incorrect arguments being supllied
- Added input validation for the params file
- Added version information to the log files

### 1.5.2
- Fixed a bug where the incorrect genes are labeled on a volcano plot
- Fixed a bug where SiSaNA required the "preprocess" parameter to be defined in the params.yml file, even if the user was not running the preprocess step.
- Added functionality for the modeProcess argument from netZooPy
- Added functionality for specifiying whether to run PANDA on CPU or GPU
- Added a hotstart option for the lioness_to_pickle_df.py step, so that if the user has an issue in the "sisana generate" step after network reconstruction, they can jump back into the analysis without needing to re-reconstruct all the networks

### 1.5.1
- Added the example Jupyter notebook for clustering
- Fixed the spacing issue and the text size of the legend in the heatmap step. Note that this fix only works when using a single metadata group (e.g. subtype). Two or more metadata types are still not supported at this time.
- Added the example output results_summarized.html file to the GitHub repo.
- Fixed the example file downloads so it downloads from the Kuijjer Lab repo 

### 1.5.0
- Added option to run LIONESS on a subset of samples for batching of the "generate" step. PANDA is still ran on all samples, so the same background is used for each batch. 
- Removed some hard-coded values in volcano plot, so now all user-defined values should work properly
- Changed the labels on the volcano plot to make it more clear which groups were which
- When visualizing heatmaps previously, an error would occur if the metadata header was not in 
the correct format. This has now been fixed.
- The "sisana compare means" command has been changed to just "sisana compare". Likewise, the "sisana compare survival" command is now just "sisana survival"
- Added a "sisana summarize" command, which summarizes the analysis the user performs in an html format.

### 1.4.1
- Changed the listed files in the log output files to all not have "./" at the beginning of their paths, just for consistency's sake
- Added more descriptive descriptions of each parameter in the params.yml file
- Added this changelog file

### 1.4.0
This new version of SiSaNA now generates log files as well, allowing the ability to find and reference the parameters you used for each analysis performed. The log files are automatically generated from the root project directory into a folder titled log_files. An example of the log file is given below.

### 1.3.0
SiSaNA now has the option to create clustermaps, allowing the user to visualize clusters of samples and genes/TFs. For this option, users also specify categorical metadata columns to color the samples on. 

### 1.2.0
With this newest version of SiSaNA, all analysis is now performed with the use of a params.yml file instead of specifying arguments directly via the command line.

### 1.1.0
Many features are now available in the new version. These include the following:

1. Filtering of all PPI, motif, and gene expression files, which is a prerequisite for running PANDA/LIONESS.
2. Filtering the output of Lioness for only edges found in the prior.
3. Calculation of the in-and out-degree of genes and TFs, respectively.
4. Reducing the number of decimal places in either the PANDA/LIONESS output or the calculated in-/out-degrees, which greatly saves on storage space.
5. Extraction of specific TFs/genes, which is useful for analyses such as limma (see Ritchie et al., 2015). An option to perform analyses with limma is not available as part of this software.
6. Comparison of groups identified in dimension reduction techniques such as UMAP or tSNE. These include comparisons of TFs/genes between two groups, survival analysis between groups, and gene set enrichment analysis (GSEA)
7. Visualization of these results via volcano plots, box plots, violin plots, and heatmaps is also possible.

### 1.0.0
Initial release. SiSaNA can reconstruct networks using PANDA/LIONESS as well as calculate in- and out-degree and compare degrees between groups.
