import pandas as pd

class NotASubsetError(Exception):
    """
    Raise when the supplied list (genes/samples) is not a subset of the items in the data frame
    
    Attributes:
        user_gene_list : List of genes/samples the user inputs
        data_genes : List of genes/samples in the dataset
        dtype : data type that makes up the list, either "genes" or "samples" 
        message -- explanation of the error
    """

    def __init__(self, user_list, data_list, dtype, message="Error: The genes in the supplied list are not a subset of the genes in the dataset."):
        
        self.user_list = user_list
        self.data_list = data_list
        self.dtype = dtype
        
        items_missing = list(set(self.user_list).difference(self.data_list)) # Find genes/samples not in the data provided
        print(items_missing)
        if self.dtype == "genes":
            print("\nError: The following genes in the supplied gene list are not present in the data provided:")
        elif self.dtype == "samples":
            print("\nError: The following samples in the supplied sample list are not present in the data provided:")
        [print(i) for i in items_missing]
        print("\n")

        self.message = message 
        super().__init__(self.message)
        
class IncorrectHeaderError(Exception):
    """
    Raise when the user's metadata file header is not in the correct format of "name,group"
    
    Attributes:
        metadf : Data frame of the user's metadata file
        message -- explanation of the error
    """
    def __init__(self, metadf=pd.DataFrame, message="Error: The header of the metadata file must be in the format 'name,group'."):
        self.metadf = metadf
        print(f"\nThe header of your metadata file looks like the following:\n {self.metadf.head()} \n")
        self.message = message 
        super().__init__(self.message)
