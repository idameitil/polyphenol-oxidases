from Bio import SeqIO
import json

fasta_sequences = SeqIO.parse('data/pfam/protein-matching-PF00264-shortheaders.fasta', 'fasta')

output_filename = 'data/pfam/protein-matching-PF00264-shortheaders-filtered.fasta'

json_list = json.load(open('data/pfam/protein-matching-PF00264.json'))

json_dict = {}
for entry in json_list:
    json_dict[entry['metadata']['accession']] = entry

with open(output_filename, 'w') as outfile:
    for fasta in fasta_sequences:
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
        if len(fasta.seq) > 100 and len(fasta.seq) < 1000 and score < 1e-5 and hit_length > 100:
            outfile.write(f">{fasta.id}\n{fasta.seq}\n")