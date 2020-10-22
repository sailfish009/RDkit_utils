import pandas as pd
from rdkit import Chem

def mol_checker(x, y):
    '''
    Takes a pandas dataframe of RDkit molecules and response variable, returns 
    a pandas dataframe with duplicate and missing values removed.
    
    Arguments:
        x: Pandas series of RDkit molecules 

        y: Pandas series of response variable

    '''

    print(f'Number of input molecules: {len(x)}')
    
    assert type(x) == pd.Series and type(y) == pd.Series, print('X and Y must be pandas series')
    
    # Identify and remove missing molecules
    missing = x.isna()
    total_missing = sum(missing)

    if total_missing > 0:
        x, y = x[-missing], y[-missing]
        print(f'Number of missing molecules removed: {total_missing}')
    else:
        print('Zero missing molecules')

    # Convert to smiles
    tmp = x.apply(lambda mol: Chem.MolToSmiles(mol)) 

    # Identify and remove duplicate molecules
    duplicates = tmp.duplicated()
    total_duplicates = sum(tmp.duplicated())
    
    if total_duplicates > 0:
        tmp, y = tmp[-duplicates], y[-duplicates]
        print(f'Number of duplicate molecules removed: {total_duplicates}')
    else:
        print('Zero duplicate molecules')
        
    x_new = tmp.apply(lambda mol: Chem.MolFromSmiles(mol)) 
    dataframe = pd.concat([x_new, y], axis=1).reset_index(drop=True, inplace=False)
    print(f'Number of output molecules: {len(dataframe)}')
    return dataframe 