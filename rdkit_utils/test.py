import pandas as pd
from transformers import start_pipeline, smiles_to_mols, mols_to_smiles
from processors import select_columns, remove_missing_samples, remove_duplicates, standardise_mols
from descriptors import calc_morgan_fp

df = pd.read_csv('../data/solubility.csv')

df.loc[1, 'smiles'] = df.loc[0, 'smiles']

clean_df = (df
.pipe(start_pipeline)
.pipe(smiles_to_mols, 'smiles')
#.pipe(standardise_mols, 'mols')
.pipe(remove_duplicates, 'mols')
.pipe(calc_morgan_fp, 'mols', y = 'solubility')
)

print(clean_df.head())
