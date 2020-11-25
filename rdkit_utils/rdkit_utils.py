import time

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import DataStructs
from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect
from rdkit.Chem.MolStandardize import Standardizer

### Wrappers ###

def log(func):
    ''' Log function for pandas pipeline functions
    '''
    def wrapper(dataframe, *args, **kwargs):
        start = time.time()
        result = func(dataframe, *args, **kwargs)
        end = time.time()
        print(f'{func.__name__} took {(end-start):.2f}s and returned {result.shape[0]} samples and {result.shape[1]} feature')
        return result
    return wrapper

def plot(func):
    ''' Plots the distribution from the results of the pipeline function
    '''
    def wrapper(dataframe, *args, **kwargs):
        results = func(dataframe, *args, **kwargs)
        sns.displot(results, kde = True)
        plt.show()
    return wrapper

### Transformers

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

@log
def standardise_mols(dataframe, column_name = 'mols', transformation = 'super_parent'):
    ''' Applies the MOL VS standardiser with a given transformation
        to the given column

        Arguments:
            dataframe: Pandas dataframe
            column_name: Column which standardiser is applied to
            transformation: Standardisation transform to be performed, see RDkit
            MOL VS for complete list of options
               '''
    std = Standardizer()
    dataframe[column_name] = dataframe[column_name].apply(lambda x: getattr(std, transformation)(x))
    return dataframe

### Selectors ###

@log
def select_columns(dataframe, *args):
    ''' Select the columns of a dataframe and all the rows

        Arguments:
            dataframe: Pandas dataframe
            *args: Column names to be selected
    '''
    return dataframe.loc[:, [*args]]

@log
def select_molw_range(dataframe, min_mw, max_mw):
    pass

@log
def select_organic(dataframe):
    pass

### Descriptors ###

@log
def calc_rdkit_descriptors():
    pass

@log
def calc_mordred_descriptors():
    pass

def morgan_fp(mol, fp_radius = 2, fp_bits = 2048):
    ''' A rdkit molecule and returns the morgan fingerprints as a np.array

    Arguments:
        mol: rdkit molecule
        fp_radius: Radius for the morgan algorithm
        fp_bits: Number of bits in the fingerprint
    '''
    fingerprint = GetMorganFingerprintAsBitVect(mol, radius = fp_radius, nBits = fp_bits)
    array = np.zeros((0,))
    DataStructs.ConvertToNumpyArray(fingerprint, array)
    return array

@log
def calc_morgan_fp(dataframe, column_name = 'mols', radius = 2, bits = 2048, y = None):
    ''' Calculates the morgan fingerprints for a column of rdkit mols

        Arguments:
            dataframe: Pandas dataframe
            column_name: Name of the column containing the rdkit mols
            radius: radius for the calculation of the fingerprints
            bits: Bit length of the fingerprints
            y: Column name of addition colum to concatenate with fingerprints
    '''
    fingerprints = [morgan_fp(mol = x, fp_radius = radius, fp_bits = bits) for x in dataframe[column_name]]
    fps = pd.DataFrame(np.vstack(fingerprints))
    if y is not None:
        return pd.concat([fps, dataframe.loc[:, y]], axis = 1)
    else:
        return fps

@log
def calc_mol2vec():
    pass

### Removers ###

@log
def remove_missing_mols(dataframe, column_name = 'mols'):
    ''' Removes rows with missing values in the given column

        Arguments:
            dataframe: Pandas dataframe
            column_name: Column which will be searhed for missing values
    '''
    missing = dataframe[column_name].isna()
    total_missing = sum(missing)
    if total_missing > 0:
        return dataframe[-missing].reset_index(drop = True, inplace = False)
    else:
        return dataframe

@log
def remove_duplicate_mols(dataframe, column_name = 'mols'):
    ''' Converts given column to smiles strings to check for duplicate molecules.
        Molecules should be standardised before removing duplciates

        Arguments:
            dataframe: Pandas dataframe
            column_name: Column which will be searched for duplicates
    '''
    if type(dataframe[column_name].values[0]) is Chem.rdchem.Mol:

        # If the type is rdkit mol convert to smiles first
        dataframe['temp_smiles'] = dataframe[column_name].apply(lambda x: Chem.MolToSmiles(x))
        duplicates = dataframe['temp_smiles'].duplicated()
        total_duplicates = sum(duplicates)

        if total_duplicates > 0:
            return dataframe[-duplicates].reset_index(drop = True, inplace = False).drop('temp_smiles', axis = 1)
        else:
            return dataframe.drop('temp_smiles', axis = 1)

    if type(dataframe[column_name].values[0]) is str:

            # If mols already as smiles strings just check for duplicates
            duplicates = dataframe[column_name].duplicated()
            total_duplicates = sum(duplicates)

            if total_duplicates > 0:
                return dataframe[-duplicates].reset_index(drop = True, inplace = False)
            else:
                return dataframe

### Visualisers ###

@plot
def plot_logp(dataframe, column_name = 'mols', ):
    ''' Calculates the average partition coefficient for each rdkit mol provided in column_name
        Arguments:
            dataframe: Pandas dataframe
            column_name: Name of the column containing the rdkit mols
    '''
    return dataframe[column_name].apply(lambda x: Descriptors.MolLogP(x))

@plot
def plot_molweight(dataframe, column_name = 'mols'):
    ''' Calculates the average Molecular weight for each rdkit mol provided in column_name

        Arguments:
            dataframe: Pandas dataframe
            column_name: Name of the column containing the rdkit mols
    '''
    return dataframe[column_name].apply(lambda x: Descriptors.MolWt(x))
