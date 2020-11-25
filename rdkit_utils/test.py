import pandas as pd
import rdkit_utils as rdu

df = pd.read_csv('../data/solubility.csv')

df_descriptors = (df
.pipe(rdu.start_pipeline)
.pipe(rdu.smiles_to_mols, 'smiles')
.pipe(rdu.remove_missing_mols)
.pipe(rdu.standardise_mols)
.pipe(rdu.remove_duplicate_mols)
.pipe(rdu.calc_morgan_fp, 'mols', y = 'solubility')
)

print(df_descriptors.head())

df_plot = (df
.pipe(rdu.start_pipeline)
.pipe(rdu.smiles_to_mols, 'smiles')
.pipe(rdu.remove_missing_mols)
.pipe(rdu.standardise_mols)
.pipe(rdu.remove_duplicate_mols)
.pipe(rdu.plot_molweight)
)
