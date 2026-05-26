import pandas as pd
import numpy as np
import sys 

class NumColorsNumGroupsMismatchError(Exception):
    """
    Raise when the user does not supply the correct number of colors for their given sub-categories (e.g. only one color was given for the category "Sex" that contained "male" and "female" for subcategories)

    Attributes:
        category: str, Name of the category 
        unique_categories: list
        unique_subcategories: int, number of given color codes for the given category 
    """
    def __init__(self, group_param_name: str, color_param_name: str, num_groups: int, num_unique_color_codes: str):
        self.group_param_name = group_param_name
        self.color_param_name = color_param_name
        self.num_groups = num_groups
        self.num_unique_color_codes = num_unique_color_codes
        self.message = f"\n\nError: The number of colors specified in your {self.color_param_name} parameter in the params.yml file ({self.num_unique_color_codes}) does not match the number of groups found in the {self.group_param_name} category ({self.num_groups}). Please fix and try running this script again.\n" 
        super().__init__(self.message)

class NotASubsetError(Exception):
    """
    Raise when the supplied list (genes/samples) is not a subset of the items in the data frame
    
    Attributes:
        user_list : List of genes/samples the user inputs
        data_list : List of genes/samples in the dataset
        dtype : data type that makes up the list, either "genes" or "samples" 
        message -- explanation of the error
    """

    def __init__(self, user_list, data_list, dtype, message="Error: The items in the supplied list are not a subset of the data."):
        
        self.user_list = user_list
        self.data_list = np.unique(data_list)
        self.dtype = dtype
                
        # Not sure if order matters here. Keeping it for testing purposes later, but it at least does catch the error for now
        # items_missing = list(set(self.user_list).difference(self.data_list)) # Find genes/samples not in the data provided
        items_missing = list(set(self.data_list).difference(self.user_list)) # Find genes/samples not in the data provided
        print(f"\nError: The following {self.dtype} in the supplied list are not present in the data provided:")
        [print(i) for i in items_missing]
        print(f"\nYou have provided the following {self.dtype}:")
        [print(i) for i in self.user_list]
        print(f"\nThe following are the {self.dtype} found in the data provided:")
        [print(i) for i in self.data_list]
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
    def __init__(self, metadf: pd.DataFrame, message="Error: The header of the metadata file must be in the format 'name,group'."):
        self.metadf = metadf
        print(f"\nThe header of your metadata file looks like the following:\n {self.metadf.head()} \n")
        self.message = message 
        super().__init__(self.message)

class WrongAmountOfColorsError(Exception):
    """
    Raise when the user does not supply the correct number of colors for their given sub-categories (e.g. only one color was given for the category "Sex" that contained "male" and "female" for subcategories)

    Attributes:
        category: str, Name of the category 
        unique_categories: list
        unique_subcategories: int, number of given color codes for the given category 
    """
    def __init__(self, category: str, num_unique_subcategories: int, num_unique_color_codes: str):
        self.category = category
        self.num_unique_subcategories = num_unique_subcategories
        self.num_unique_color_codes = num_unique_color_codes
        self.message = f"\n\nError: The number of colors specified in your category_column_colors parameter in the params.yml file ({self.num_unique_color_codes}) does not match the number of unique subcategories found in your metadata file for the {category} category ({self.num_unique_subcategories}). Please fix and try running this script again.\n" 
        super().__init__(self.message)

class TooManyCoresError(Exception):
    """
    Raise when the user supplies more cores than they have samples for in the dataset for the "ncores" parameter in the params.yml file for the reconstruct command)

    Attributes:
        category: str, Name of the category 
        unique_categories: list
        unique_subcategories: int, number of given color codes for the given category 
    """
    def __init__(self, ncores: str, nsamps: int):
        self.ncores = ncores
        self.nsamps = nsamps
        self.message = f"\n\nError: You have requested more cores ({self.ncores}) than you have samples ({self.nsamps}). Please ensure 'ncores' <= number of samples.\n"
        super().__init__(self.message)
        
class ExtensionMismatchError(Exception):
    """
    Raise when the user supplies a file with an extension that does not match the delimiter they specified in the params.yml file.

    Attributes:
        category: str, Name of the category 
    """
    def __init__(self, filename: str, ext: str):
        self.filename = filename
        self.ext = ext
        self.message = f"\n\nError: The supplied file ({filename}) has a different extension than the one indicated in the params file ({ext}). Please ensure the data format matches the extension you indicated in the params file and that you have set the parameters for it correctly.\n"
        super().__init__(self.message)
        
