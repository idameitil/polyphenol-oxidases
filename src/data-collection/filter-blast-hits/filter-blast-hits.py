import pandas as pd
from Bio import SeqIO
import sys

evalue_threshold = float(sys.argv[1])

unique_hits_filename = 'data/blast/unique-hits.tsv'

df = pd.read_csv(unique_hits_filename, sep='\t', names=['accession', 'score', 'evalue'])

df_evaluefilter = df[df.evalue < evalue_threshold]

with open(f'data/blast/unique-hits-{evalue_threshold}.fasta', 'w') as outfile:
    with open(f'data/blast/unique-hits-{evalue_threshold}-length150-1000.fasta', 'w') as outfile2:
        fasta_sequences = SeqIO.parse('data/blast/unique-hits.fasta', 'fasta')
        for fasta in fasta_sequences:
            id, seq = fasta.id, str(fasta.seq)
            if id in list(df_evaluefilter.accession):
                outfile.write(f">{id}\n{seq}\n")
                if len(seq) > 150 and len(seq) < 1000:
                    outfile2.write(f">{id}\n{seq}\n")