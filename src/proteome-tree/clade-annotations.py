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

seeds2species = {
    'AbPPO1': 'Agaricus_bisporus',
    'AbPPO2': 'Agaricus_bisporus',
    'AbPPO3': 'Agaricus_bisporus',
    'AbPPO4': 'Agaricus_bisporus',
    'AbPPO5': 'Agaricus_bisporus',
    'AbPPO6': 'Agaricus_bisporus',
    'AoMelB': 'Aspergillus_oryzae',
    'AoCO4': 'Aspergillus_oryzae',
    'HsTrp1': 'Homo_sapiens',
    'MdPPO1': 'Malus_domestica',
}

uniprot_seeds = {
    "P17643": "HsTrp1",
    "A0A1S9DJX3": "AoCO4",
    "C7FF04": "AbPPO3",
    "C7FF05": "AbPPO4",
    "A0A8H7C4Z8": "AbPPO6",
    "A0A8H7C551": "AbPPO1",
    "A0A8H7F0X6": "AbPPO3",
    "A0A8H7F114": "AbPPO2",
    "A0A8H7F118": "AbPPO4",
    "A0A8H7F178": "AbPPO5",
    "A0A1S9DK56": "AoMelB",
    "A0A498IXV7": "MdPPO1",
    "D5DZK6": "BmTyr"
}

seeds = [
'ScTyr',
'NpsF',
'GriF',
'BtTyr',
'SaTyr',
'VsTyr',
'BmTyr',
'VvPPO',
'SlPPO1',
'MdPPO1',
'JrPPO1',
'CgAUS1',
'IbCO',
'EdHc',
'RvHc',
'HdHc',
'MgHc',
'TpHc',
'CgPPO-473',
'CgPPO-266',
'MtPPO-809',
'TtPPO',
'PpPPO-c2092',
'MtPPO7',
'HjTyr',
'PsTyr'
]

def get_clade(acc):
    if acc in acc2clade:
        return acc2clade[acc]
    else:
        return 's'

def add_to_dict(species, clade, species2clades):
    if species not in species2clades:
        species2clades[species] = {clade: 1}
    else:
        if clade not in species2clades[species]:
            species2clades[species][clade] = 1
        else:
            species2clades[species][clade] += 1
    return species2clades

def make_df(included_accessions, acc2clade):
    accs_done = []
    species2clades = {}
    for acc in included_accessions:
        if acc == "D5DZK6":
            continue
        elif acc in seeds2species:
            clade = get_clade(acc)
            species = seeds2species[acc]
        elif acc in seeds:
            continue
        else:
            # exclude if the protein has been there before (proteins with several domains)
            acc_stripped = acc.split('/')[0]
            if acc_stripped in accs_done:
                continue
            accs_done.append(acc_stripped)
            # get species
            mask = df_uniprot_hits.protein_accession == acc_stripped
            species = df_uniprot_hits[mask].species.item()
            species = species.replace(' ', '_')
            # get clade
            clade = get_clade(acc)
        species2clades = add_to_dict(species, clade, species2clades)
    df = pd.DataFrame(species2clades).fillna(0).T
    df.index.name = 'species'
    return df

# Run with all
included_accessions = get_included_accessions('data/epa-ng/filtered-out/ref-query-linsi.fa')
acc2clade = get_clade_members(clade_dir='data/epa-ng/filtered-out/clades')
df = make_df(included_accessions, acc2clade)
df.to_csv('data/epa-ng/filtered-out/clades.csv')

# Run with all
included_accessions = get_included_accessions('data/proteome-tree/sequences-fungi-with-lignin-degraders-order-notFiltered.trimmed.fa') + ['AoCO4', 'AoMelB', 'AbPPO1', 'AbPPO2', 'AbPPO3', 'AbPPO4', 'AbPPO5', 'AbPPO6']
acc2clade = get_clade_members(clade_dir='data/epa-ng/fungi-with-lignin-degraders/clades')
df = make_df(included_accessions, acc2clade)
df.to_csv('data/epa-ng/fungi-with-lignin-degraders/clades.csv')
