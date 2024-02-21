import pandas as pd
import math

df = pd.read_csv('data/eggnog/MM_7t_157_w.emapper.annotations.tsv', sep='\t', skiprows=[0,1,2,3])

data = {}
max = 11
for index, row in df.iterrows():
    acc = row['#query']
    if not pd.isna(row['eggNOG_OGs']):
        OGs = row['eggNOG_OGs'].split(',')
        data[acc] = OGs + [''] * (max - len(OGs))

df2 = pd.DataFrame(data).transpose()

df2.to_csv('data/eggnog/OGs.tsv', sep='\t')