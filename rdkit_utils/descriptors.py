from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect
import pandas as pd
import numpy as np
from utils import log
from utils import morgan_fp

### Descriptors ###

@log
def calc_morgan_fp(dataframe, column_name, radius = 2, bits = 2048, y = None):
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
        
### Fingerprints ###


### Mol2Vec ###



