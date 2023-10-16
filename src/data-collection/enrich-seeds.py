import pandas as pd
from common import get_taxon

input_filename = "data/seeds.tsv"

df_seed_table = pd.read_csv(input_filename, sep='\t', encoding='latin-1')

desired_ranks = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
for rank in desired_ranks:
    df_seed_table[rank] = df_seed_table.taxon.apply(lambda x: get_taxon(x, rank))

# Get interproscan data
interproscan_output_filename = 'data/seeds.interproscan'
df_interproscan = pd.read_csv(interproscan_output_filename, sep='\t', names=['protein_accession', 'md5', 'sequence_length', \
                                                                'analysis', 'signature_accession', 'signature_discription',\
                                                                'start_location', 'stop_location', 'score', 'status',\
                                                                'date', 'interpro_annotations_accession', \
                                                                'interpro_annotations_description'])
pd.set_option('display.max_columns', None)
unique_pfams = df_interproscan[df_interproscan.analysis == 'Pfam'].signature_accession.unique()
data = {}
for index, row in df_seed_table.iterrows():
    data[row.id] = []
    print(row.id)
    for pfam_name in unique_pfams:
        hits = df_interproscan[(df_interproscan.protein_accession == row.id) & (df_interproscan.signature_accession == pfam_name)]
        if hits.empty:
            data[row.id].append('')
        else:
            data[row.id].append(f"{hits.start_location.values[0]}-{hits.stop_location.values[0]}")
df_pfam_positions = pd.DataFrame.from_dict(data, orient='index', columns=unique_pfams)

# Merge
df_merged = pd.merge(left=df_seed_table, right=df_pfam_positions, how='left', left_on='id', right_index=True)

# Write output file
output_filename = "data/seeds-enriched.tsv"
df_merged.to_csv(output_filename, sep='\t', index=False)