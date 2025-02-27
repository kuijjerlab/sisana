from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import pandas as pd
import numpy as np
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt
import seaborn as sns
import os 

def scale_data(df):
    '''
    Function that scales data to be used for PCA decomposition later

    Parameters 
    ----------
        df : panda data frame
            A df of just numeric values
    '''   
    scaler = StandardScaler()
    scaled = scaler.fit_transform(df)
    return(scaled)


def pca_fit_transform(scaled, nc, df_index):
    '''
    Function that performs a PCA decomposition on scaled data and returns a dictionary of the calculated components

    Parameters 
    ----------
        scaled : panda data frame
            A df of scaled values (output of scale_data())
        nc : int
            The number of components to keep from the PCA decomposition. Must match the same number used in pca_cum_var()
        df_index : list
            A list of sample names to give the output data frame as index
    '''
    
    pca = PCA(n_components = nc)
    pca_table = pca.fit_transform(scaled)
    # print("pca table:")
    # print(pca_table)
    colnames = [f"PC{x}" for x in range(1, nc+1)]
    pca_df = pd.DataFrame(data = pca_table, columns = colnames)
    pca_df = pca_df.set_index(df_index)
    pca_dict = {}
    pca_dict["df"] = pca_df
    pca_dict["pcs"] = pca_table
    return(pca_dict)

def plot_pca(df, catcol, output_path):
    '''
    Function that visualizes the resulting PCA plot on the first two PCs

    Parameters 
    ----------
        df : panda data frame
            A data frame including the metadata (sample categories) and PCs for each sample
        catcol : str
            Choice of column to use for assigning categories in plot
        output_path : str
            Path to file name to save the figure as
    '''
    # plt.figure(figsize=(5,5))
    # print(df)
    # contin_vals = df[continuous_col]
    # print(contin_vals)
    
    # cmap = sns.cubehelix_palette(rot=-.2, as_cmap=True)
    
    plot = sns.scatterplot(x="PC1", y="PC2",
                        hue = catcol,
                        palette=sns.color_palette(),
                        data=df,
                        legend="full",
                        alpha=0.5) 


    fig = plot.get_figure()
    
    fig.savefig(output_path, bbox_inches='tight')
    fig.show()
    plt.clf()

def km(pca_df, kmax, numsamps, clusttype):
    '''
    Function that performs kmeans clustering on data for multiple values of k and returns the silhouette scores for each
    
    Adapted from https://realpython.com/k-means-clustering-python/

    Parameters 
    ----------
        pca_df : panda data frame
            A df of PCA values output by pca_fit_transform
        kmax : int
            The maximum value of k, where k is the number of clusters
        numsamps : int
            The number of samples in the data frame
        clusttype : str [choices are kmeans, hkmeans, or hierarchical]
            The type of clustering to perform. kmeans performs standard kmeans clustering while hkmeans uses the Hartigan method 
    ''' 
    
    silhouette_scores = []
    kmax_new = kmax + 1

    if clusttype == "kmeans":
        kmeans_kwargs = {
            "init": "random",
            "n_init": 10,
            "max_iter": 300,
            "random_state": 42,
        }

        for k in range(2, kmax_new):
            kmeans = KMeans(n_clusters = k, **kmeans_kwargs)
            kmeans.fit(pca_df)
            
            if len(np.unique(kmeans.labels_)) <= numsamps - 1: 
                score = silhouette_score(pca_df, kmeans.labels_)
                silhouette_scores.append(score)
                
    elif clusttype == "hkmeans":
        from sklearn.feature_extraction.text import TfidfVectorizer
        from sklearn import metrics
        from hkmeans import HKMeans
        
        for k in range(2, kmax_new):
            hkmeans = HKMeans(n_clusters = k, random_state = 128, n_init = 10,
                            n_jobs = 16, max_iter = 15)
            hkmeans.fit(pca_df)
            
            if len(np.unique(hkmeans.labels_)) <= numsamps - 1: 
                # print(hkmeans.labels_)
                score = silhouette_score(pca_df, hkmeans.labels_)
                silhouette_scores.append(score)
    
    elif clusttype == "hierarchical":
        from scipy.cluster.hierarchy import ward, fcluster
        
        for k in range(2, kmax_new):
            w = ward(pca_df) # Perform Ward's linkage      
            fl = fcluster(w, k, criterion='maxclust')
                
            if len(np.unique(fl)) <= numsamps - 1: 
                score = silhouette_score(pca_df, fl)
                silhouette_scores.append(score)
         
        # if len(np.unique(fl)) <= numsamps - 1: 
        #         score = silhouette_score(pca_df, fl)
        #         silhouette_scores.append(score)
    
    return(silhouette_scores)

def assign_to_cluster(pca_df, data_t, k, clusttype):
    '''
    Function that assigns the samples to clusters and returns a data frame containing the cluster assignment for each sample
    
    Parameters 
    ----------
        pca_df : panda data frame
            A df of PCA values that were calculated by pca_fit_transform
        data_t : panda data frame
            A transposed data frame, where samples are rows and variables (e.g. genes) are columns.
        k : int
            The number of clusters 
        clusttype : str [choices are kmeans, hkmeans, hierarchical, or consensus]
            The type of clustering to perform. kmeans performs standard kmeans clustering while hkmeans uses the Hartigan method 
    '''     
    cluster_map = pd.DataFrame()
    cluster_map['sample'] = data_t.index.values

    if clusttype == "kmeans":
        kmeans_kwargs = {
            "init": "random",
            "n_init": 10,
            "max_iter": 300,
            "random_state": 42,
        }
        
        km = KMeans(n_clusters=k, **kmeans_kwargs).fit(pca_df)
        cluster_map['cluster'] = km.labels_

    elif clusttype == "hkmeans":
        from sklearn.feature_extraction.text import TfidfVectorizer
        from sklearn import metrics
        from hkmeans import HKMeans
        
        hkmeans = HKMeans(n_clusters = k, random_state = 128, n_init = 10,
                        n_jobs = 16, max_iter = 15)
        hkmeans.fit(pca_df)
        cluster_map['cluster'] = hkmeans.labels_
        
    elif clusttype == "hierarchical":
        from scipy.cluster.hierarchy import ward, fcluster
        
        w = ward(pca_df) # Perform Ward's linkage      
        fl = fcluster(w, k, criterion='maxclust')
        cluster_map['cluster'] = fl
        
    elif clusttype == "consensus":
        from .consensusClustering import ConsensusCluster
        
        cc = ConsensusCluster(KMeans, 2, k, 500)
        cc.fit(pca_df)
        # consensus_res = cc.predict()
        res = cc.predict_data(pca_df)
        score = silhouette_score(pca_df, res)
        cluster_map['cluster'] = res

        # print(cluster_map)
        # print(score)
        
    print("\nBelow, the first column is the cluster number and the second column is the number of samples assigned to that cluster.")
    print(cluster_map['cluster'].value_counts())
    print("\n")
    return(cluster_map)
