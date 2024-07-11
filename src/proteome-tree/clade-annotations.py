import os
import pandas as pd
from Bio import SeqIO

df_uniprot_hits = pd.read_csv('data/pfam/protein-matching-PF00264-interproscan2.tsv', sep='\t', low_memory=False)

def read_file(filename):
    with open(filename) as file:
        lines = [line.rstrip() for line in file]
    return lines

def get_included_accessions(fasta_filename):
    fasta_sequences = SeqIO.parse(fasta_filename, 'fasta')
    return [fasta.id for fasta in fasta_sequences]

def get_clade_members(clade_dir):
    acc2clade = {}
    for clade in os.listdir(clade_dir):
        members = read_file(f"{clade_dir}/{clade}") 
        for acc in members:
            acc2clade[acc] = clade
    return acc2clade

def make_df(included_accessions, acc2clade):
    accs_done = []
    species2clades = {}
    count = 0
    for acc in included_accessions:
        count += 1
        # exclude seeds
        if '.' not in acc:
            continue
        if '/' in acc:
            acc = acc.replace('/', '0')
        acc_stripped = acc.split('.')[0]
        # exclude if the protein has been there before (proteins with several domains)
        if acc_stripped in accs_done:
            continue
        accs_done.append(acc_stripped)
        mask = df_uniprot_hits.protein_accession == acc_stripped
        if len(df_uniprot_hits[mask]) == 0:
            continue
        species = df_uniprot_hits[mask].species.item()
        species_underscore = species.replace(' ', '_')
        if acc in acc2clade:
            clade = acc2clade[acc]
        else:
            clade = 'singletons'
        clade_short = clade[0]
        if species_underscore not in species2clades:
            species2clades[species_underscore] = {clade_short: 1}
        else:
            if clade_short not in species2clades[species_underscore]:
                species2clades[species_underscore][clade_short] = 1
            else:
                species2clades[species_underscore][clade_short] += 1
    df = pd.DataFrame(species2clades).fillna(0).T
    df.index.name = 'species'
    return df

# Run with all
# included_accessions = get_included_accessions('data/proteome-tree/all-one_proteome_per_class.fa')
included_accessions = get_included_accessions('data/epa-ng/filtered-out/ref-query-linsi.fa')
acc2clade = get_clade_members(clade_dir='data/mrbayes/all-seeds-0619/clades/members')
df = make_df(included_accessions, acc2clade)
df.to_csv('data/mrbayes/all-seeds-0619/clades/clades.csv')

# Run with fungi one per order
included_accessions = get_included_accessions(fasta_filename='data/epa-ng/fungi-order/query.fa')
acc2clade = get_clade_members(clade_dir='data/fungi-genome-figure/clades/members')
df = make_df(included_accessions, acc2clade)
df.to_csv('data/fungi-genome-figure/clades.csv')
