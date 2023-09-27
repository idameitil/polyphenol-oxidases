import pandas as pd
import random

# def make_domain_color_file(df):
#     outfilename = f"data/itol-label-files/domain.txt"
#     domain2color = {'bacterial': '#880808', 'plant': '#00FFFF', 'fungal': '#008000', 'animal': '#FFFF00'}
#     with open(outfilename, "w") as file:
#         header = f"DATASET_COLORSTRIP\nSEPARATOR COMMA\nDATASET_LABEL,domain\nCOLOR,#ff0000\nDATA\n"
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

def make_activity_label_file(df):
    outfilename = f"data/itol-label-files/activity.txt"
    activity2color = {'tyrosinase': '#880808', 'catechol oxidase': '#00FFFF', 'Hemocyanin': '#008000', 'Pro-phenoloxidase': '#FFFF00', 'O-aminophenol oxidase': '#00FF00'}
    with open(outfilename, "w") as file:
        header = f"DATASET_COLORSTRIP\nSEPARATOR COMMA\nDATASET_LABEL,activity\nCOLOR,#ff0000\nDATA\n"
        file.write(header)
        for index, row in df.iterrows():
            if row.enzyme_class not in activity2color:
                continue
            file.write(f"{row.id},{activity2color[row.enzyme_class]},{row.enzyme_class}\n")

def make_methoxylation_label_file(df):
    outfilename = f"data/itol-label-files/methoxylation.txt"
    methoxylation2color = {0: '#FF0000', 1: '#00FFFF', 2: '#00FF00'}
    with open(outfilename, "w") as file:
        header = f"DATASET_COLORSTRIP\nSEPARATOR COMMA\nDATASET_LABEL,methoxylation\nCOLOR,#ff0000\nDATA\n"
        file.write(header)
        for index, row in df.iterrows():
            if row.methoxylation not in methoxylation2color:
                continue
            file.write(f"{row.id},{methoxylation2color[row.methoxylation]},{row.methoxylation}\n")

def make_yes_no_label_file(df, column_name):
    outfilename = f"data/itol-label-files/{column_name}.txt"
    value2color = {'Yes': '#008000', 'No': '#000000'}
    with open(outfilename, "w") as file:
        header = f"DATASET_COLORSTRIP\nSEPARATOR COMMA\nDATASET_LABEL,{column_name}\nCOLOR,#ff0000\nDATA\n"
        file.write(header)
        for index, row in df.iterrows():
            if row[column_name] not in value2color:
                continue
            file.write(f"{row.id},{value2color[row[column_name]]},{row[column_name]}\n")

seed_df = pd.read_csv('data/seeds-enriched.tsv', sep='\t')
# make_domain_color_file(seed_df)
make_taxonomy_label_file(seed_df)
make_activity_label_file(seed_df)
# make_methoxylation_label_file(seed_df)
make_yes_no_label_file(seed_df, '0_methoxylations')
make_yes_no_label_file(seed_df, '1_methoxylation')
make_yes_no_label_file(seed_df, '2_methoxylations')
make_yes_no_label_file(seed_df, 'Short_fungal')
make_yes_no_label_file(seed_df, 'Long_fungal')
make_yes_no_label_file(seed_df, 'Monophenolase_activity')
make_yes_no_label_file(seed_df, 'Diphenolase_activity')
make_yes_no_label_file(seed_df, 'Nitrosation_activity')
make_yes_no_label_file(seed_df, 'Tioether bond')