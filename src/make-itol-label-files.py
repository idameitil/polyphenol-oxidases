import pandas as pd
import random

# def make_domain_color_file(df):
#     outfilename = f"data/itol-label-files/domain.txt"
#     domain2color = {'bacterial': '#880808', 'plant': '#00FFFF', 'fungal': '#008000', 'animal': '#FFFF00'}
#     with open(outfilename, "w") as file:
#         header = f"DATASET_COLORSTRIP\nSEPARATOR COMMA\nDATASET_LABEL,pubmed\nCOLOR,#ff0000\nDATA\n"
#         file.write(header)
#         for index, row in df.iterrows():
#             file.write(f"{row.id},{domain2color[row.domain]},{row.domain}\n")

def make_taxonomy_label_file(df):
    wanted_ranks = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    for rank in wanted_ranks:
        outfilename = f"data/itol-label-files/{rank}.txt"
        with open(outfilename, "w") as file:
            header = f"DATASET_COLORSTRIP\nSEPARATOR TAB\nDATASET_LABEL\t{rank}\nCOLOR\t#ff0000\nDATA\n"
            file.write(header)
            tax2color = dict()
            for index, row in df.iterrows():
                acc, tax = row['id'], row[rank]
                if tax not in tax2color:
                    color = '#' + "%06x" % random.randint(0, 0xFFFFFF)
                    tax2color[tax] = color
                file.write(f"{acc}\t{tax2color[tax]}\t{tax}\n")

seed_df = pd.read_csv('data/seeds-enriched.tsv', sep='\t')
# make_domain_color_file(seed_df)
make_taxonomy_label_file(seed_df)