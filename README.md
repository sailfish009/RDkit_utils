# RDkit_utils

Utility script for the transformation, processing and visualisation of molecular data for cheminformatics 
and machine learning. 

# Installation

Clone the directory using git or download the single python file:  

```
git clone https://github.com/prbradshaw/RDkit_utils.git

-- or --

wget https://raw.githubusercontent.com/prbradshaw/RDkit_utils/main/rdkit_utils/rdkit_utils.py

```
If cloning the whole repo, create a virtual environment using conda

```
conda env create --file venv.txt --name venv

```

## Example
The package relies on pandas pipe functionality which allow the chaining of functions on a pandas dataframe 
or series.The example below loads in a csv file containing smiles strings and converts them to rdkit mols. 
The mols are then standardised, duplicates are removed and the morgan fingerprint for each molecule is calculated.

```python3
import pandas as pd
import rdkit_utils as rdu

df = pd.read_csv('...')

clean_df = (df
.pipe(rdu.start_pipeline)
.pipe(rdu.smiles_to_mols, 'smiles')
.pipe(rdu.standardise_mols, 'mols')
.pipe(rdu.remove_duplicates, 'mols')
.pipe(rdu.calc_morgan_fp, 'mols', y = 'solubility')
)

```

Work in progress...
