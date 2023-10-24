import pandas as pd
import random

def make_taxonomy_label_files(df, blast_hits=False):
    wanted_ranks = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    for rank in wanted_ranks:
        # Get output filename
        if blast_hits:
            output_filename = f"data/itol-label-files/{rank}-blast-hits.txt"
        else:
            output_filename = f"data/itol-label-files/{rank}-seeds.txt"
        # Write file
        with open(output_filename, "w") as file:
            header = f"DATASET_COLORSTRIP\nSEPARATOR TAB\nDATASET_LABEL\t{rank}\nCOLOR\t#ff0000\nDATA\n"
            file.write(header)
            tax2color = dict()
            for index, row in df.iterrows():
                acc, tax = row['protein_accession'], row[rank]
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
            file.write(f"{row.protein_accession},{activity2color[row.enzyme_class]},{row.enzyme_class}\n")

def make_methoxylation_label_file(df):
    outfilename = f"data/itol-label-files/methoxylation.txt"
    methoxylation2color = {0: '#FF0000', 1: '#00FFFF', 2: '#00FF00'}
    with open(outfilename, "w") as file:
        header = f"DATASET_COLORSTRIP\nSEPARATOR COMMA\nDATASET_LABEL,methoxylation\nCOLOR,#ff0000\nDATA\n"
        file.write(header)
        for index, row in df.iterrows():
            if row.methoxylation not in methoxylation2color:
                continue
            file.write(f"{row.protein_accession},{methoxylation2color[row.methoxylation]},{row.methoxylation}\n")

def make_yes_no_label_file(df, column_name):
    outfilename = f"data/itol-label-files/{column_name}.txt"
    value2color = {'Yes': '#008000', 'No': '#000000'}
    with open(outfilename, "w") as file:
        header = f"DATASET_COLORSTRIP\nSEPARATOR COMMA\nDATASET_LABEL,{column_name}\nCOLOR,#ff0000\nDATA\n"
        file.write(header)
        for index, row in df.iterrows():
            if row[column_name] not in value2color:
                continue
            file.write(f"{row.protein_accession},{value2color[row[column_name]]},{row[column_name]}\n")

def make_aguilera_subclass_label_file(df):
    outfilename = f"data/itol-label-files/Aguilera_subclass.txt"
    value2color = {'alpha': '#FF0000', 'beta': '#00FFFF', 'gamma': '#00FF00'}
    with open(outfilename, "w") as file:
        header = f"DATASET_COLORSTRIP\nSEPARATOR COMMA\nDATASET_LABEL,Aguilera_subclass\nCOLOR,#ff0000\nDATA\n"
        file.write(header)
        for index, row in df.iterrows():
            if row['Aguilera_subclass'] not in value2color:
                continue
            file.write(f"{row.protein_accession},{value2color[row['Aguilera_subclass']]},{row['Aguilera_subclass']}\n")

def make_aguilera_subclass_label_file_text(df):
    outfilename = f"data/itol-label-files/Aguilera_subclass_text.txt"
    value2color = {'alpha': '#FF0000', 'beta': '#00FFFF', 'gamma': '#00FF00'}
    with open(outfilename, "w") as file:
        header = f"DATASET_TEXT\nSEPARATOR COMMA\nDATASET_LABEL,Aguilera_subclass\nCOLOR,#ff0000\nDATA\n"
        file.write(header)
        for index, row in df.iterrows():
            if row['Aguilera_subclass'] not in value2color:
                continue
            file.write(f"{row.protein_accession},{row['Aguilera_subclass']},-1,{value2color[row['Aguilera_subclass']]},bold,1,0\n")

def make_domain_label_file(df, blast_hits=False):
    if blast_hits:
        output_filename = "data/itol-label-files/domain-blast-hits.txt"
    else:
        output_filename = "data/itol-label-files/domain-seeds.txt"
    with open(output_filename, "w") as file:
        header = f"DATASET_DOMAINS\nSEPARATOR COMMA\nDATASET_LABEL,Domains\nCOLOR,#ff0000\nDATA\n"
        file.write(header)
        domain2color = dict()
        domains = [column_name for column_name in df.columns if column_name.startswith('domain_')]
        for index, row in df.iterrows():
            try:
                file.write(f"{row.protein_accession},{len(row.sequence)}")
            except:
                continue
            for domain_name in domains:
                domain_description = domain_name.split('_')[2]
                if domain_name not in domain2color:
                    color = '#' + "%06x" % random.randint(0, 0xFFFFFF)
                    domain2color[domain_name] = color
                if not pd.isnull(row[domain_name]):
                    # Several occurences of this domain in the sequence
                    if ',' in row[domain_name]:
                        occurences = row[domain_name].split(',')
                        for occurence in occurences:
                            start, stop = occurence.split('-')
                            file.write(f",OC|{start}|{stop}|{domain2color[domain_name]}|{domain_description}")
                    # One occurence of the domain in the sequence
                    else:
                        start, stop = row[domain_name].split('-')
                        file.write(f",OC|{start}|{stop}|{domain2color[domain_name]}|{domain_description}")
            file.write('\n')

# Make seed label files
seed_df = pd.read_csv('data/seeds-enriched.tsv', sep='\t')
# make_domain_color_file(seed_df)
make_taxonomy_label_files(seed_df)
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
make_domain_label_file(seed_df)

# Make blast hit label files
df_blast_hits = pd.read_csv('data/blast/unique-hits-1e-60-length150-1000-cd-hit65-enriched.tsv', sep='\t')
make_domain_label_file(df_blast_hits, blast_hits=True)
make_taxonomy_label_files(df_blast_hits, blast_hits=True)