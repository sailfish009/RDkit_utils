import pandas as pd

# Dataset from https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/delaney-processed.csv

df = pd.read_csv('delaney-processed.csv')
# Select smiles strings and response variable only
df = df.loc[:, ['measured log solubility in mols per litre', 'smiles']]
df = df.rename(columns = {'measured log solubility in mols per litre': 'solubility'})
df.to_csv('solubility.csv', index = False)

