import pandas as pd
from Bio import SeqIO

# Seeds
interproscan_output_filename = 'data/seeds.interproscan'
fasta_filename = 'data/seeds.fa'

fasta2seq = {}
fasta_sequences = SeqIO.parse(open(fasta_filename), 'fasta')
for fasta in fasta_sequences:
    fasta2seq[fasta.id] = str(fasta.seq)

df = pd.read_csv(interproscan_output_filename, sep='\t', names=['protein_accession', 'md5', 'sequence_length', \
                                                                'analysis', 'signature_accession', 'signature_discription',\
                                                                'start_location', 'stop_location', 'score', 'status',\
                                                                'date', 'interpro_annotations_accession', \
                                                                'interpro_annotations_description'])

df_PPO_domains = df[(df['signature_accession']=='PF00264')|(df['signature_accession']=='PF00372')]

output_filename = 'data/seeds-pfam-domains.fa'
accessions_done = []
with open(output_filename, 'w') as outfile:
    for index, row in df_PPO_domains.iterrows():
        if row.protein_accession in accessions_done:
            continue
        pfam_domain = fasta2seq[row.protein_accession][row.start_location-1:row.stop_location-1]
        outfile.write(f">{row.protein_accession}\n{pfam_domain}\n")
        accessions_done.append(row.protein_accession)

# Blast hits
# enriched_hits_df = pd.read_csv("data/blast/unique-hits-1e-60-length150-1000-cd-hit65-enriched.tsv", sep='\t')
# print(enriched_hits_df[~enriched_hits_df["domain_PF00264_Common central domain of tyrosinase"].isna()\
#       |~enriched_hits_df["domain_PF00372_Hemocyanin, copper containing domain"].isna()])

# output_filename = 'data/blast/unique-hits-1e-60-length150-1000-cd-hit65-pfam-domains.fa'
# with open(output_filename, 'w') as outfile:
#     for index, row in enriched_hits_df.iterrows():
#         if not pd.isnull(row["domain_PF00264_Common central domain of tyrosinase"]):
#             if ',' in row["domain_PF00264_Common central domain of tyrosinase"]:
#                 continue
#             else:
#                 start, stop = row["domain_PF00264_Common central domain of tyrosinase"].split('-')
#                 outfile.write(f">{row.protein_accession}\n{}")