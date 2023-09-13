import pandas as pd

def make_domain_color_file(df):
    outfilename = f"data/itol-label-files/domain.txt"
    domain2color = {'bacterial': '#880808', 'plant': '#00FFFF', 'fungal': '#008000', 'animal': '#FFFF00'}
    with open(outfilename, "w") as file:
        header = f"DATASET_COLORSTRIP\nSEPARATOR COMMA\nDATASET_LABEL,pubmed\nCOLOR,#ff0000\nDATA\n"
        file.write(header)
        for index, row in df.iterrows():
            file.write(f"{row.id},{domain2color[row.domain]},{row.domain}\n")

seed_df = pd.read_csv('data/seeds.tsv', sep='\t')
make_domain_color_file(seed_df)