import pandas as pd
from common import get_taxon, read_interproscan_output, make_interproscan_df

# Read seed table
input_filename = "data/seeds.tsv"
df_seed_table = pd.read_csv(input_filename, sep='\t', encoding='latin-1')

# Get taxonomy
desired_ranks = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
for rank in desired_ranks:
    df_seed_table[rank] = df_seed_table.taxon.apply(lambda x: get_taxon(x, rank))

# Get interproscan data
filename = f"data/seeds.interproscan"
acc2domain, unique_domains = read_interproscan_output(filename)
df_interproscan = make_interproscan_df(acc2domain, unique_domains)

# Merge
merged_df = pd.merge(df_seed_table, df_interproscan, how='left', on='protein_accession')
merged_df.to_csv('data/seeds-enriched.tsv', sep='\t', index=False)