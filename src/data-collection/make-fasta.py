import pandas as pd

df = pd.read_csv('data/seeds.tsv', sep='\t', encoding='latin-1')

output_filename = 'data/seeds.fa'
with open(output_filename, 'w') as outfile:
    for index, row in df.iterrows():
        outfile.write(f">{row.protein_accession}\n{row.sequence}\n")
