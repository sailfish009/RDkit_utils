import matplotlib.pyplot as plt
import seaborn as sns

from rdkit.Chem import Descriptors

def plot_wrapper(func):
    ''' Plots a histogram from the results of the pipeline function
    '''
    def wrapper(dataframe, *args, **kwargs):
        results = func(dataframe, *args, **kwargs)
        sns.distplot(results)
        plt.show()
    return wrapper

@plot_wrapper
def log_p(dataframe, column_name):
    ''' Calculates the average partition coefficient for each rdkit mol provided in column_name 
        
        Arguments:
            dataframe: Pandas dataframe
            column_name: Name of the column containing the rdkit mols
    '''
    return dataframe[column_name].apply(lambda x: Descriptors.MolLogP(x))

@plot_wrapper
def mol_weight(dataframe, column_name):
    ''' Calculates the average Molecular weight for each rdkit mol provided in column_name 
        
        Arguments:
            dataframe: Pandas dataframe
            column_name: Name of the column containing the rdkit mols
    '''
    return dataframe[column_name].apply(lambda x: Descriptors.MolWt(x))

