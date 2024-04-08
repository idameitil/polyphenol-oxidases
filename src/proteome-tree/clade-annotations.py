import os
import pandas as pd
clade_dir = 'data/proteome-tree/clades'

def read_file(filename):
    with open(filename) as file:
        lines = [line.rstrip() for line in file]
    return lines

clade2members = {}
accs = []
acc2clade = {}
for clade in os.listdir(clade_dir):
    members = read_file(f"{clade_dir}/{clade}") 
    clade2members[clade] = members
    accs.extend(members)
    for acc in members:
        acc2clade[acc] = clade

df_uniprot_hits = pd.read_csv('data/pfam/protein-matching-PF00264-interproscan2.tsv', sep='\t')

species2clades = {}
for acc in accs:
    mask = df_uniprot_hits.protein_accession == acc
    if len(df_uniprot_hits[mask]) == 0:
        continue
    species = df_uniprot_hits[mask].species.item()
    clade = acc2clade[acc]
    if species not in species2clades:
        species2clades[species] = {clade: 1}
    else:
        if clade not in species2clades[species]:
            species2clades[species][clade] = 1
        else:
            species2clades[species][clade] += 1

df = pd.DataFrame(species2clades).fillna(0).T
df.index.name = 'species'

df.to_csv('data/proteome-tree/clades/clades.csv')
