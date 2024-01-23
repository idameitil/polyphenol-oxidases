from Bio import SeqIO
from common import get_taxon

fasta_filename = 'data/pfam/protein-matching-PF00264.fasta'

fasta_sequences = SeqIO.parse(fasta_filename, 'fasta')

output_filename = "data/pfam/protein-matching-PF00264-fungi.fasta"
with open(output_filename, 'w') as outfile:
    for entry in fasta_sequences:
        taxid = entry.description.split(':')[-1]
        kingdom = get_taxon(taxid, 'kingdom')
        if kingdom == "Fungi":
            outfile.write(f">{entry.description}\n{entry.seq}\n")