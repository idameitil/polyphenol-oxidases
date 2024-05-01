import os
import pandas as pd
from Bio import SeqIO
clade_dir = 'data/mrbayes/all/clades/members'

def read_file(filename):
    with open(filename) as file:
        lines = [line.rstrip() for line in file]
    return lines

def get_included_accessions():
    filename = 'data/proteome-tree/all-one_proteome_per_class.fa'
    fasta_sequences = SeqIO.parse(filename, 'fasta')
    return [fasta.id for fasta in fasta_sequences]

included_accessions = get_included_accessions()

clade2members = {}
acc2clade = {}
for clade in os.listdir(clade_dir):
    members = read_file(f"{clade_dir}/{clade}") 
    clade2members[clade] = members
    for acc in members:
        acc2clade[acc] = clade

df_uniprot_hits = pd.read_csv('data/pfam/protein-matching-PF00264-interproscan2.tsv', sep='\t', low_memory=False)

species2clades = {}
for acc in included_accessions:
    mask = df_uniprot_hits.protein_accession == acc
    if len(df_uniprot_hits[mask]) == 0:
        continue
    species = df_uniprot_hits[mask].species.item()
    if acc in acc2clade:
        clade = acc2clade[acc]
    else:
        clade = 'singletons'
    if species not in species2clades:
        species2clades[species] = {clade: 1}
    else:
        if clade not in species2clades[species]:
            species2clades[species][clade] = 1
        else:
            species2clades[species][clade] += 1

df = pd.DataFrame(species2clades).fillna(0).T
df.index.name = 'species'

df.to_csv('data/mrbayes/all/clades/clades.csv')
