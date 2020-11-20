# RDkit_utils

Package for the transformation, processing and visualisation of molecular data for cheminformatics and machine learning. 

The package is based on RDkit and pandas and allows:

- Removal of 
- Removal of duplicate molecules
- Standardisation of molecular graphs
- Calculation of molecular descriptors and fingerprints


## Example
The package relies on pandas pipe functionality which allow the chaining of functions on a pandas dataframe or series.The 
example below loads in a csv file, converts the smiles strings to rdkit mol, standardises the mols, removes any duplicates
and calculates the morgan fingerprints for each molecule.  

```python3
import pandas as pd

rom transformers import start_pipeline, smiles_to_mols
from processors import remove_duplicates, standardise_mols
from descriptors import calc_morgan_fp

df = pd.read_csv('...')

clean_df = (df
.pipe(start_pipeline)
.pipe(smiles_to_mols, 'smiles')
.pipe(standardise_mols, 'mols')
.pipe(remove_duplicates, 'mols')
.pipe(calc_morgan_fp, 'mols', y = 'solubility')
)

```

Working in progress...
