import pandas as pd
from rdkit import Chem
from utils import log
 
@log
def start_pipeline(dataframe):
    ''' Starts a pandas pipeline by copying the dataframe to avoid overriding
        the original dataframe.

        Arguments:
            dataframe: Pandas dataframe
    '''
    return dataframe.copy()

@log
def smiles_to_mols(dataframe, smiles_column_name, rename_column = 'mols'):
    ''' Takes a columns of smiles strings and transforms them to rdkit mols 
        Arguments:
            dataframe: Pandas dataframe
            smiles_column_name: The column of the dataframe containing the smile strings
            rename_column: Rename column, if original column is desired set to False
    '''
    dataframe[smiles_column_name] = dataframe[smiles_column_name].apply(lambda x: Chem.MolFromSmiles(x))
    
    if rename_column is not False:
        return dataframe.rename(columns = {smiles_column_name: rename_column})
    else:
        return dataframe
    
@log
def mols_to_smiles(dataframe, mol_column_name, rename_column = 'smiles'):
    ''' Takes a columns of rdkit mols and returns a column of smiles strings  
        Arguments:
            dataframe: Pandas dataframe
            mol_column_name: The column of the dataframe containing the rdkit mol
            rename_column: Rename column, if original column is desired set to False
    '''
    dataframe[mol_column_name] = dataframe[mol_column_name].apply(lambda x: Chem.MolToSmiles(x))
    
    if rename_column is not False:
        return dataframe.rename(columns = {mol_column_name: rename_column})
    else:
        return dataframe
