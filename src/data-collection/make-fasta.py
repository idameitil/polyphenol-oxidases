import pandas as pd

df = pd.read_csv('data/seeds.tsv', sep='\t')

output_filename = 'data/seeds.fasta'
with open(output_filename, 'w') as outfile:
    for index, row in df.iterrows():
        outfile.write(f">{row.id}\n{row.sequence}\n")