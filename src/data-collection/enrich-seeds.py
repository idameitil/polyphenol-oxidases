import pandas as pd
from common import get_taxon

input_filename = "data/seeds.tsv"
output_filename = "data/seeds-enriched.tsv"

df = pd.read_csv(input_filename, sep='\t')

desired_ranks = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
for rank in desired_ranks:
    df[rank] = df.taxon.apply(lambda x: get_taxon(x, rank))

df.to_csv(output_filename, sep='\t', index=False)