from rdkit import Chem
from rdkit.Chem.MolStandardize import Standardizer
from utils import log

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

### Removers ###

@log
def remove_missing_samples(dataframe, column_name):
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
def remove_duplicates(dataframe, column_name):
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

### Standardiser ###

@log
def standardise_mols(dataframe, column_name, transformation = 'super_parent'):
    ''' Applies the MOL VS standardiser with a given transformation
        to the given column

        Arguments:
            dataframe: Pandas dataframe
            transformation: Standardisation transform to be performed, see RDkit
            MOL VS for complete list of options
            column_name: Column which standardiser is applied to, Default = 'ROMol'
    '''
    std = Standardizer()
    dataframe[column_name].apply(lambda x: getattr(std, transformation)(x))
    return dataframe


