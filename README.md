# RDkit_utils

Package for the transformation, processing and visualisation of molecular data for cheminformatics and machine learning. 

The package is based on RDkit and pandas, features include:

- Removal of missing molecules
- Removal of duplicate molecules
- Standardisation of molecular graphs
- Calculation of molecular descriptors and fingerprints
- Plotting distrubution of molecular descriptors and fingerprints

## Example
The package relies on pandas pipe functionality which allow the chaining of functions on a pandas dataframe or series.The 
example below loads in a csv file containing smiles strings and converts them to rdkit mols.  The mols are then standardised, duplicates 
are removed and the morgan fingerprint for each molecule is calculated.

```python3
import pandas as pd

from rdkit_utils.transformers import start_pipeline, smiles_to_mols
from rdkit_utils.processors import remove_duplicates, standardise_mols
from rdkit_utils.descriptors import calc_morgan_fp

df = pd.read_csv('...')

clean_df = (df
.pipe(start_pipeline)
.pipe(smiles_to_mols, 'smiles')
.pipe(standardise_mols, 'mols')
.pipe(remove_duplicates, 'mols')
.pipe(calc_morgan_fp, 'mols', y = 'solubility')
)

```

Work in progress...
