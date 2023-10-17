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
                if rank == 'kingdom':
                    kingdom2color = {'Viridiplantae': '#00FF00', 'nan': '#FF0000', 'Fungi': '#964B00', 'Metazoa': '#ffff00'}
                    if tax not in kingdom2color:
                        color = '#FF0000'
                    else:
                        color = kingdom2color[tax]
                    file.write(f"{acc}\t{color}\t{tax}\n")
                else:
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

def make_aguilera_subclass_label_file(df):
    outfilename = f"data/itol-label-files/Aguilera_subclass.txt"
    value2color = {'alpha': '#FF0000', 'beta': '#00FFFF', 'gamma': '#00FF00'}
    with open(outfilename, "w") as file:
        header = f"DATASET_COLORSTRIP\nSEPARATOR COMMA\nDATASET_LABEL,Aguilera_subclass\nCOLOR,#ff0000\nDATA\n"
        file.write(header)
        for index, row in df.iterrows():
            if row['Aguilera_subclass'] not in value2color:
                continue
            file.write(f"{row.id},{value2color[row['Aguilera_subclass']]},{row['Aguilera_subclass']}\n")

def make_aguilera_subclass_label_file_text(df):
    outfilename = f"data/itol-label-files/Aguilera_subclass_text.txt"
    value2color = {'alpha': '#FF0000', 'beta': '#00FFFF', 'gamma': '#00FF00'}
    with open(outfilename, "w") as file:
        header = f"DATASET_TEXT\nSEPARATOR COMMA\nDATASET_LABEL,Aguilera_subclass\nCOLOR,#ff0000\nDATA\n"
        file.write(header)
        for index, row in df.iterrows():
            if row['Aguilera_subclass'] not in value2color:
                continue
            file.write(f"{row.id},{row['Aguilera_subclass']},-1,{value2color[row['Aguilera_subclass']]},bold,1,0\n")

def make_pfam_label_file(df):
    outfilename = f"data/itol-label-files/pfam.txt"
    pfam_names = ['PF00264', 'PF14830', 'PF03723', 'PF03722', 'PF00372', 'PF12143', 'PF12142', 'PF18132']
    pfam2color = {'PF00264':'#8BCED0', 'PF14830':'#8D0083', 'PF03723':'#F7F700', 'PF03722':'#B7C8B8', \
                  'PF00372':'#002E91', 'PF12143':'#E73841', 'PF12142':'#FFD9D9', 'PF18132':'#F1B40D'}
    with open(outfilename, "w") as file:
        header = f"DATASET_DOMAINS\nSEPARATOR COMMA\nDATASET_LABEL,Domains\nCOLOR,#ff0000\nDATA\n"
        file.write(header)
        for index, row in df.iterrows():
            file.write(f"{row.id},{len(row.sequence)}")
            for pfam_name in pfam_names:
                if not pd.isnull(row[pfam_name]):
                    # Several occurences of this pfam in the sequence
                    if ',' in row[pfam_name]:
                        occurences = row[pfam_name].split(',')
                        for occurence in occurences:
                            start, stop = occurence.split('-')
                            file.write(f",OC|{start}|{stop}|{pfam2color[pfam_name]}|{pfam_name}")
                    # One occurence of the pfam in the sequence
                    else:
                        start, stop = row[pfam_name].split('-')
                        file.write(f",OC|{start}|{stop}|{pfam2color[pfam_name]}|{pfam_name}")
            file.write('\n')

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
make_aguilera_subclass_label_file(seed_df)
make_aguilera_subclass_label_file_text(seed_df)
make_pfam_label_file(seed_df)