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

def make_binary_label_files(df):
    outfilename = f"data/itol-label-files/binary.txt"
    shape2number = {'rectangle': 1, 'circle': 2, 'star':3, 'right_pointing_triangle':4, \
                    'left_pointing_triangle':5, 'check_mark': 6}
    columns = {'0_methoxylations': ['rectangle', '#00FFFF'], 
               '1_methoxylation': ['rectangle', '#FF0000'], 
               '2_methoxylations': ['rectangle', '#008000'], 
               'Short_fungal': ['left_pointing_triangle', '#00FFFF'], 
               'Long_fungal': ['left_pointing_triangle', '#FF0000'],
               'Monophenolase_activity': ['circle', '#00FFFF'],
               'Diphenolase_activity': ['circle', '#FF0000'], 
               'Nitrosation_activity': ['star', '#000000'], 
               'Tioether bond': ['check_mark', '#000000']}
    with open(outfilename, "w") as file:
        # Write header
        header = f"DATASET_BINARY\nSEPARATOR COMMA\nDATASET_LABEL,binary\nCOLOR,#ff0000\n"
        file.write(header)
        # Write field shapes
        file.write("FIELD_SHAPES")
        for column in columns:
            shape = columns[column][0]
            file.write(f",{shape2number[shape]}")
        file.write('\n')
        # Write field labels
        file.write("FIELD_LABELS")
        for column in columns:
            file.write(f",{column}")
        file.write('\n')        
        # Write field colours
        file.write("FIELD_COLORS")
        for column in columns:
            color = columns[column][1]
            file.write(f",{color}")
        file.write('\n')
        # Write data
        file.write('DATA\n')
        for index, row in df.iterrows():
            file.write(f"{row.protein_accession}")
            for column_name in columns:
                if row[column_name] == 'Yes':
                    fill = 1
                elif row[column_name] == 'No':
                    fill = 0
                else:
                    fill = -1
                file.write(f",{fill}")
            file.write('\n')

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
            if row['Subclass'] not in value2color:
                continue
            if pd.isnull(row.protein_accession):
                continue
            file.write(f"{row.protein_accession},{row['Subclass']},-1,{value2color[row['Subclass']]},bold,1,0\n")

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
                if blast_hits:
                    file.write(f"{row.protein_accession},{len(row.seq)}")
                else:
                    file.write(f"{row.protein_accession},{len(row.sequence)}")
            except:
                continue
            for domain_name in domains:
                domain_description = '_'.join(domain_name.split('_')[1:3])
                if domain_description == 'NON_CYTOPLASMIC' or domain_description == 'CYTOPLASMIC_DOMAIN' \
                    or domain_description == 'SIGNAL_PEPTIDE':
                    continue
                if domain_description.startswith('PF'):
                    shape = 'HH'
                elif domain_description.startswith('SignalP'):
                    shape = 'RE'
                else:
                    shape = 'EL'
                if domain_name not in domain2color:
                    color = '#' + "%06x" % random.randint(0, 0xFFFFFF)
                    domain2color[domain_name] = color
                if not pd.isnull(row[domain_name]):
                    # Several occurences of this domain in the sequence
                    if ',' in row[domain_name]:
                        occurences = row[domain_name].split(',')
                        for occurence in occurences:
                            start, stop = occurence.split('-')
                            file.write(f",{shape}|{start}|{stop}|{domain2color[domain_name]}|{domain_description}")
                    # One occurence of the domain in the sequence
                    else:
                        start, stop = row[domain_name].split('-')
                        file.write(f",{shape}|{start}|{stop}|{domain2color[domain_name]}|{domain_description}")
            file.write('\n')

# Make seed label files
seed_df = pd.read_csv('data/seeds-enriched.tsv', sep='\t')
make_taxonomy_label_files(seed_df)
make_activity_label_file(seed_df)
make_binary_label_files(seed_df)
make_aguilera_subclass_label_file(seed_df)
# make_aguilera_subclass_label_file_text(seed_df)
make_domain_label_file(seed_df)

# Make blast hit label files
df_blast_hits = pd.read_csv('data/blast/unique-hits-1e-60-length150-1000-cd-hit65-enriched.tsv', sep='\t')
make_domain_label_file(df_blast_hits, blast_hits=True)
make_taxonomy_label_files(df_blast_hits, blast_hits=True)

# Make aguilera label files
df_aguilera = pd.read_excel('data/Aguilera-data/aguilera_with_seq.xlsx')
make_aguilera_subclass_label_file_text(df_aguilera)