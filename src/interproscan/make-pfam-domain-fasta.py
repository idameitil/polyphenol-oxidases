import pandas as pd
from Bio import SeqIO

interproscan_output_filename = 'data/seeds.interproscan'
fasta_filename = 'data/seeds.fa'

fasta2seq = {}
fasta_sequences = SeqIO.parse(open(fasta_filename), 'fasta')
for fasta in fasta_sequences:
    fasta2seq[fasta.id] = str(fasta.seq)

df = pd.read_csv(interproscan_output_filename, sep='\t', names=['protein_accession', 'md5', 'sequence_length', \
                                                                'analysis', 'signtature_accession', 'signature_discription',\
                                                                'start_location', 'stop_location', 'score', 'status',\
                                                                'date', 'interpro_annotations_accession', \
                                                                'interpro_annotations_description'])

df_PPO_domains = df[(df['signtature_accession']=='PF00264')|(df['signtature_accession']=='PF00372')]

output_filename = 'data/seeds-pfam-domains.fa'
accessions_done = []
with open(output_filename, 'w') as outfile:
    for index, row in df_PPO_domains.iterrows():
        if row.protein_accession in accessions_done:
            continue
        pfam_domain = fasta2seq[row.protein_accession][row.start_location-1:row.stop_location-1]
        outfile.write(f">{row.protein_accession}\n{pfam_domain}\n")
        accessions_done.append(row.protein_accession)