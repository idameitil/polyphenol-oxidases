from Bio import SeqIO
import json

fasta_sequences = SeqIO.parse('data/pfam/protein-matching-PF00264.fasta', 'fasta')

output_filename = 'data/pfam/protein-matching-PF00264-filtered20.fasta'

json_list = json.load(open('data/pfam/protein-matching-PF00264.json'))

threshold = 1e-25

json_dict = {}
for entry in json_list:
    json_dict[entry['metadata']['accession']] = entry
count = 0
with open(output_filename, 'w') as outfile:
    for fasta in fasta_sequences:
        accession = fasta.id.split('|')[0]
        if len(json_dict[accession]['entries']) > 1:
            print(f"{accession} has several entries")
            continue
        for entry in json_dict[accession]['entries'][0]['entry_protein_locations']:
            if len(entry['fragments']) > 1:
                print(f"{accession} has several fragments")
                continue
            score = entry['score']
            start = entry['fragments'][0]['start']
            end = entry['fragments'][0]['end']
            if score < threshold:
                outfile.write(f">{accession}/{start}-{end}\n{fasta.seq}\n")
            else:
                count += 1
print(count)