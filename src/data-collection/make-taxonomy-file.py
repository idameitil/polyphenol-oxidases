from Bio import SeqIO
from common import get_taxon
import pandas as pd

fasta_filename = 'data/pfam/protein-matching-PF00264.fasta'

fasta_sequences = SeqIO.parse(fasta_filename, 'fasta')

acc2info = {}
for entry in fasta_sequences:
    taxid = entry.description.split(':')[-1]
    kingdom = get_taxon(taxid, 'kingdom')
    acc = entry.id.split('|')[0]
    acc2info[acc] = {'taxid': taxid, 'description': entry.description}
    # Get taxonomy
    desired_ranks = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    for rank in desired_ranks:
        acc2info[acc][rank] = get_taxon(taxid, rank)
    acc2info[acc]['seq'] = entry.seq

df = pd.DataFrame.from_dict(acc2info, orient='index')
df = df.reset_index().rename(columns={"index":"protein_accession"})
df.to_csv("data/pfam/protein-matching-PF00264.tsv", sep='\t')