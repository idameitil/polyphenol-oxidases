import pandas as pd

df = pd.read_csv('data/seeds.tsv', sep='\t', encoding='latin-1')

output_filename = 'data/seeds.fa'
with open(output_filename, 'w') as outfile:
    for index, row in df.iterrows():
        outfile.write(f">{row.descriptive_name}\n{row.full_sequence}\n")

# df = pd.read_excel('data/Aguilera-data/aguilera_with_seq.xlsx')
# output_filename = 'data/Aguilera-data/aguilera.fa'
# with open(output_filename, 'w') as outfile:
#     for index, row in df.iterrows():
#         if pd.isnull(row.seq):
#             continue
#         outfile.write(f">{row.protein_accession}\n{row.seq}\n")
