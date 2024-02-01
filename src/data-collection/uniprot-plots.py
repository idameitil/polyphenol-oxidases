from Bio import SeqIO
import matplotlib.pyplot as plt
import json
import numpy as np

fasta_sequences = SeqIO.parse('data/pfam/protein-matching-PF00264-shortheaders.fasta', 'fasta')

json_list = json.load(open('data/pfam/protein-matching-PF00264.json'))

json_dict = {}
for entry in json_list:
    json_dict[entry['metadata']['accession']] = entry

lenghts = []
coverages = []
scores = []
for fasta in fasta_sequences:
    lenghts.append(len(fasta.seq))

    if len(json_dict[fasta.id]['entries']) > 1:
        print(f"{fasta.id} has several entries")
        continue
    if len(json_dict[fasta.id]['entries'][0]['entry_protein_locations']) > 1:
        print(f"{fasta.id} has several entry_protein_locations")
        continue
    if len(json_dict[fasta.id]['entries'][0]['entry_protein_locations'][0]['fragments']) > 1:
        print(f"{fasta.id} has several fragments")
        continue
    coverages.append(json_dict[fasta.id]['entries'][0]['entry_protein_locations'][0]['fragments'][0]['end'] - json_dict[fasta.id]['entries'][0]['entry_protein_locations'][0]['fragments'][0]['start'])
    scores.append(json_dict[fasta.id]['entries'][0]['entry_protein_locations'][0]['score'])

plt.hist(lenghts, bins=50, range=[0,1000])
plt.xlabel('Length')
plt.savefig('data/pfam/plots/length.png')
plt.clf()

plt.hist(coverages, bins=50)
plt.xlabel('Coverage')
plt.savefig('data/pfam/plots/coverage.png')
plt.clf()

log_scores = np.log(scores)
plt.hist(log_scores, bins=100)
plt.xlabel('log(score)')
plt.savefig('data/pfam/plots/score.png')