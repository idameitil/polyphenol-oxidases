from Bio import SeqIO
import json
import pandas as pd
import sys

#########################################################
# Select proteomes
#########################################################
df = pd.read_csv('data/proteome-tree/proteome-data.tsv', sep='\t')

fungal_orders = df[df['kingdom'] == 'Fungi']['order'].unique()
selected_proteome_taxids = []
for order in fungal_orders:
    if pd.isna(order):
        continue
    best_proteome = df[df.order == order]['score'].idxmax()
    best_proteome_taxid = df.iloc[[best_proteome]].taxid.item()
    selected_proteome_taxids.append(best_proteome_taxid)

#########################################################
# Get accessions for sequences in selected proteomes
#########################################################
json_list = json.load(open('data/pfam/protein-matching-PF00264.json'))

sequences_from_selected_proteomes = list()
for entry in json_list:
    if int(entry['metadata']['source_organism']['taxId']) in selected_proteome_taxids:
        sequences_from_selected_proteomes.append(entry['metadata']['accession'])

#########################################################
# Filter and write fasta
#########################################################
json_dict = {}
for entry in json_list:
    json_dict[entry['metadata']['accession']] = entry

fasta_sequences = SeqIO.parse('data/pfam/protein-matching-PF00264-shortheaders.fasta', 'fasta')
output_filename = 'data/proteome-tree/fungal-one_proteome_per_order.fa'
with open(output_filename, 'w') as outfile:
    for fasta in fasta_sequences:
        if fasta.id not in sequences_from_selected_proteomes:
            continue
        if len(json_dict[fasta.id]['entries']) > 1:
            print(f"{fasta.id} has several entries")
            continue
        if len(json_dict[fasta.id]['entries'][0]['entry_protein_locations']) > 1:
            print(f"{fasta.id} has several entry_protein_locations")
            continue
        if len(json_dict[fasta.id]['entries'][0]['entry_protein_locations'][0]['fragments']) > 1:
            print(f"{fasta.id} has several fragments")
            continue
        hit_length = json_dict[fasta.id]['entries'][0]['entry_protein_locations'][0]['fragments'][0]['end'] - json_dict[fasta.id]['entries'][0]['entry_protein_locations'][0]['fragments'][0]['start']
        score = json_dict[fasta.id]['entries'][0]['entry_protein_locations'][0]['score']
        if len(fasta.seq) > 150 and len(fasta.seq) < 1000 and score < 1e-5 and hit_length > 150:
            outfile.write(f">{fasta.id}\n{fasta.seq}\n")
