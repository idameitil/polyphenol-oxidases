import pandas as pd
from common import read_interproscan_output, make_interproscan_df

# Read enriched unique hits
df_unique_hits_enriched = pd.read_csv("data/pfam/protein-matching-PF00264.tsv", sep = '\t')

# Get interproscan data
N = 20
unique_domains = set()
acc2domain = {}
for i in range(0, N):
    filename = f"data/interproscan-uniprot/chunk{'%02d' % i}.interproscan"
    acc2domain, unique_domains = read_interproscan_output(filename, acc2domain, unique_domains)
df_interproscan = make_interproscan_df(acc2domain, unique_domains)

# Merge
merged_df = pd.merge(df_unique_hits_enriched, df_interproscan, how='left', on='protein_accession')
merged_df.to_csv('data/pfam/protein-matching-PF00264-interproscan.tsv', sep='\t', index=False)