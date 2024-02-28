import json
from common import get_taxon
from Bio import SeqIO

json_list = json.load(open('data/pfam/protein-matching-PF00264.json'))

json_dict = {}
for entry in json_list:
    json_dict[entry['metadata']['accession']] = entry

#########################################################
# Make Orthofinder fastas
#########################################################
# Get sequences from selected fungal species
genera_proteomes = list()
sequences_from_selected_fungal_species = list()
sequence2genustaxid = dict()
selected_fungal_species = dict()
for entry in json_list:
    species = entry['metadata']['source_organism']['scientificName']
    genus = entry['metadata']['source_organism']['scientificName'].split(' ')[0]
    if entry['extra_fields']['counters']['proteome'] == 1:
        if get_taxon(entry['metadata']['source_organism']['taxId'], 'kingdom') == 'Fungi':
            if genus not in genera_proteomes:
                selected_fungal_species[genus] = species
            if selected_fungal_species[genus] == species:
                sequences_from_selected_fungal_species.append(entry['metadata']['accession'])
                # sequence2order[entry['metadata']['accession']] = order
                sequence2genustaxid[entry['metadata']['accession']] = f"{genus}_{entry['metadata']['source_organism']['taxId']}"
        genera_proteomes.append(genus)
        
# Directory for Orthofinder
outdir = 'data/orthofinder'
fasta_sequences = SeqIO.parse('data/pfam/protein-matching-PF00264-shortheaders.fasta', 'fasta')
for fasta in fasta_sequences:
    if fasta.id not in sequences_from_selected_fungal_species:
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
    # if len(fasta.seq) > 100 and len(fasta.seq) < 1000 and score < 1e-5 and hit_length > 100:
    if len(fasta.seq) > 150 and len(fasta.seq) < 1000 and score < 1e-25 and hit_length > 150:
        genustaxid = sequence2genustaxid[fasta.id]
        output_filename = f"{outdir}/{genustaxid}.fa"
        with open(output_filename, 'a') as outfile:
            outfile.write(f">{fasta.id}\n{fasta.seq}\n")