import time
import numpy as np
from rdkit.Chem import DataStructs
from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect

def log(func):
    ''' Log function for pandas pipeline functions
    '''
    def wrapper(dataframe, *args, **kwargs):
        start = time.time()
        result = func(dataframe, *args, **kwargs)
        end = time.time()
        print(f'{func.__name__} took {(end - start):.2f}s and returned {result.shape[0]} samples and {result.shape[1]} features')
        return result
    return wrapper
    
def morgan_fp(mol, fp_radius = 2, fp_bits = 2048):
    ''' Takes a list or pandas series of RDkit molecule and 
    returns the morgan fingerprints as an 
    
    Arguments:
        mol: An RDkit molecules 
        radius: Radius for the morgani algorithm, Default = 2 
        bits: Number of bits in the fingerprint
    '''
    fingerprint = GetMorganFingerprintAsBitVect(mol, radius = fp_radius, nBits = fp_bits)
    array = np.zeros((0,))
    DataStructs.ConvertToNumpyArray(fingerprint, array)
    return array

