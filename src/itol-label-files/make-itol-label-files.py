import pandas as pd
import random
import json
import math
from Bio import SeqIO
from colour import Color

random.seed(10)
outdir = "data/itol-label-files"

def make_value2color(values):
    value2color = {}
    for value in values:
        if value not in value2color:
            color = '#' + "%06x" % random.randint(0, 0xFFFFFF)
            value2color[value] = color
    return value2color

def write_colour_strip_file(outfile_name, label, ids, values, value2color={}):
    if value2color == {}:
        value2color = make_value2color(values)
    with open(outfile_name, 'w') as file:
        header = f"DATASET_COLORSTRIP\nSEPARATOR TAB\nDATASET_LABEL\t{label}-strip\nCOLOR\t#ff0000\nDATA\n"
        file.write(header)
        for id,value in zip(ids,values):
            if value not in value2color:
                color = '#FFFFFF'
            else:
                color = value2color[value]
            file.write(f"{id}\t{color}\t{value}\n")

def write_colour_text_file(outfile_name, label, ids, values, value2color={}):
    if value2color == {}:
        value2color = make_value2color(values)
    with open(outfile_name, 'w') as file:
        header = header = f"DATASET_TEXT\nSEPARATOR COMMA\nDATASET_LABEL,{label}-text\nCOLOR,#000000\nDATA\n"
        file.write(header)
        for id,value in zip(ids,values):
            if value not in value2color:
                color = '#000000'
            else:
                color = value2color[value]
            file.write(f"{id},{value},-1,{color},bold,1,0\n")

def write_heatmap_text_file(outfile_name, label, df, id_name, field_names):
    with open(outfile_name, 'w') as file:
        header = f"DATASET_HEATMAP\nSEPARATOR SPACE\nDATASET_LABEL {label}-heatmap\nCOLOR #A020F0\nFIELD_LABELS {' '.join(field_names)}\nDATA\n"
        file.write(header)
        print(df)
        for index, row in df.iterrows():
            string = ''
            for field_name in field_names:
                string += ' ' + str(row[field_name])
            file.write(f"{row[id_name].replace(' ', '_')}{string}\n")

df = pd.read_csv('data.csv')
write_heatmap_text_file(f'{outdir}/clade-heatmap.txt', 'clade', df, 'species', ['a_chordata','b_plants','c_cnidaria', 'd_long_fungal','e_mollusc','f_oomycota','g_cnidaria2','h_zoopago1','i_short_fungal','j_zoopago2','k_bacteria','undefined'])

def make_taxonomy_label_files(df):
    wanted_ranks = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    for rank in wanted_ranks:
        ids = df.protein_accession.tolist()
        values = df[rank].tolist()
        if rank == 'kingdom':
            value2color = {'Viridiplantae': '#00FF00', 'nan': '#FFFFFF', 'Fungi': '#964B00', 'Metazoa': '#ffff00'}
        else:
            value2color = {}
        output_filename_text = f"data/itol-label-files/uniprot-{rank}-text.txt"
        write_colour_text_file(output_filename_text, rank, ids, values, value2color)
        output_filename_strip = f"data/itol-label-files/uniprot-{rank}-strip.txt"
        write_colour_strip_file(output_filename_strip, rank, ids, values, value2color)

def make_taxonomy_files_species_tree(df):
    wanted_ranks = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus']
    for rank in wanted_ranks:
        ids = df['species'].str.replace(' ', '_').tolist()
        values = df[rank].tolist()
        if rank == 'kingdom':
            value2color = {'Viridiplantae': '#00FF00', 'nan': '#FFFFFF', 'Fungi': '#964B00', 'Metazoa': '#ffff00'}
        else:
            value2color = {}
        output_filename_text = f"data/itol-label-files/species-tree-{rank}-text.txt"
        write_colour_text_file(output_filename_text, rank, ids, values, value2color)
        output_filename_strip = f"data/itol-label-files/species-tree-{rank}-strip.txt"
        write_colour_strip_file(output_filename_strip, rank, ids, values, value2color)

def get_included_accessions(fasta_filename):
    return [fasta.id for fasta in SeqIO.parse(fasta_filename, 'fasta')]

def make_value2first_id(ids, values, included_accessions):
    value2first_id = {}
    for id,value in zip(ids,values):
        if id not in included_accessions:
            continue
        if value not in value2first_id:
            value2first_id[value] = id
    return value2first_id

def write_arrow_file(outfile_name, label, ids, values, included_accessions, value2color={}):
    if value2color == {}:
        value2color = make_value2color(values)
    value2first_id = make_value2first_id(ids, values, included_accessions)
    with open(outfile_name, 'w') as file:
        header = header = f"DATASET_CONNECTION\nSEPARATOR COMMA\nDATASET_LABEL,{label}-arrow\nCOLOR,#00FFFF\nDATA\n"
        file.write(header)
        for id,value in zip(ids,values):
            if id not in included_accessions:
                continue
            if value not in value2color:
                color = '#000000'
            else:
                color = value2color[value]
            if value2first_id[value] == id:
                continue
            file.write(f"{value2first_id[value]},{id},2,{color},dashed,{value}\n")

def make_taxonomy_arrow_files(df, domain):
    wanted_ranks = ['class']
    for rank in wanted_ranks:
        included_accessions = get_included_accessions(f"data/proteome-tree/{domain}-one_proteome_per_{rank}.fa")
        output_filename = f"data/itol-label-files/{domain}-{rank}-arrows.txt"
        ids = df['protein_accession'].tolist()
        values = df[rank].tolist()
        write_arrow_file(output_filename, rank, ids, values, included_accessions)

def make_activity_label_file(df):
    ids = df['protein_accession'].tolist()
    values = df['enzyme_class']
    outfilename = f"data/itol-label-files/activity.txt"
    value2color = {'tyrosinase': '#880808', 'catechol oxidase': '#00FFFF', 'Hemocyanin': '#008000', 'Pro-phenoloxidase': '#FFFF00', 'O-aminophenol oxidase': '#00FF00'}
    write_colour_text_file(outfilename, 'activity', ids, values, value2color)
    
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

def make_aguilera_subclass_files(df):
    outfilename = f"data/itol-label-files/Aguilera_subclass.txt"
    value2color = {'alpha': '#FF0000', 'beta': '#00FFFF', 'gamma': '#00FF00'}
    ids = df['protein_accession'].to_list()
    values = df['Aguilera_subclass'].to_list()
    write_colour_text_file(outfilename, 'aguilera_subclass', ids, values, value2color)
    write_colour_strip_file(outfilename, 'aguilera_subclass', ids, values, value2color)

def write_spectrum_strip_file(outfile_name, label, ids, values, value2color={}, min='', max ='', max_cutoff = ''):
    if min == '' and max == '':
        min = min(values)
        max = max(values)
    n = int(max-min)
    colors = list(Color("black").range_to(Color("white"),n))
    with open(outfile_name, "w") as file:
        header = f"DATASET_COLORSTRIP\nSEPARATOR TAB\nDATASET_LABEL\t{label}\nCOLOR\t#FFC0CB\nDATA\n"
        file.write(header)
        for id, value in zip(ids, values):
            if max_cutoff != '':
                if value > max_cutoff:
                    continue
            if value > max:
                color = colors[math.floor(max)].hex
            elif value < min:
                color = colors[math.floor(min)].hex
            else:
                color = colors[math.floor(value)].hex
            file.write(f"{id}\t{color}\t{value}\n")

def make_score_label_file():
    json_list = json.load(open('data/pfam/protein-matching-PF00264.json'))
    ids, values = [], []
    for entry in json_list:
        ids.append(entry['metadata']['accession'])
        values.append(math.log10(entry[['entries'][0]['entry_protein_locations'][0]['score']]))
    output_filename = "data/itol-label-files/HMM_score.txt"
    write_spectrum_strip_file(output_filename, 'score', ids, values, min=-140, max=-5)

def make_coverage_label_file():
    json_list = json.load(open('data/pfam/protein-matching-PF00264.json'))
    ids, values = [], []
    for entry in json_list:
        ids.append(entry['metadata']['accession'])
        values.append(entry[['entries'][0]['entry_protein_locations'][0]['fragments'][0]['end'] - json_dict[acc]['entries'][0]['entry_protein_locations'][0]['fragments'][0]['start']])
    output_filename = "data/itol-label-files/coverage175.txt"
    write_spectrum_strip_file(output_filename, 'coverage', ids, values, min=16, max=405)

def make_number_of_copies_file():
    df = pd.read_csv('data/proteome-tree/proteome-data.tsv', sep='\t')
    ids, values = [], []
    for index, row in df.iterrows():
        try:
            species = row.species.replace(' ', '_')
        except:
            continue
        ids.append(species)
        values.append(row.count_tyrosinases)
    values = df.count_tyrosinases.to_list()
    output_filename = "data/itol-label-files/copies.txt"
    write_spectrum_strip_file(output_filename, 'copies', ids, values, min=0, max=30)

def make_match_length_file():
    fasta_filename = 'data/pfam/PF00264.alignment.uniprot-cleaned-filtered-withoutgaps.fa'
    entries = SeqIO.parse(fasta_filename, 'fasta')
    max = 214
    # min_match_length = 42
    min = 150
    ids,values = [], []
    for entry in entries:
        ids.appen(entry.id.split('.')[0])
        match_length = len(entry.seq)
        values.append(match_length-min)
    output_filename = "data/itol-label-files/match-length180.txt"
    write_spectrum_strip_file(output_filename, 'match_length180', ids, values, min=0, max=30, max_cutoff=180)

def make_domain_label_file(df, blast_hits=False, uniprot_hits=False):
    if blast_hits:
        output_filename = "data/itol-label-files/domain-blast-hits.txt"
    elif uniprot_hits:
        output_filename = "data/itol-label-files/domain-uniprot-hits.txt"
    else:
        output_filename = "data/itol-label-files/domain-seeds.txt"
    with open(output_filename, "w") as file:
        header = f"DATASET_DOMAINS\nSEPARATOR COMMA\nDATASET_LABEL,Domains\nCOLOR,#ff0000\nDATA\n"
        file.write(header)
        domain2color = dict()
        domains = [column_name for column_name in df.columns if column_name.startswith('domain_')]
        for index, row in df.iterrows():
            try:
                if blast_hits or uniprot_hits:
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
                elif domain_description.startswith('TMhelix'):
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

def make_OG_files():
    wanted_ranks = [0,11]
    df = pd.read_csv('data/eggnog/OGs.tsv', sep='\t', header=0, index_col=0)
    for rank in wanted_ranks:
        output_filename = f"data/itol-label-files/OG_{rank}.txt"
        # Write file
        with open(output_filename, "w") as file:
            header = f"DATASET_COLORSTRIP\nSEPARATOR TAB\nDATASET_LABEL\tOG_{rank}\nCOLOR\t#ff0000\nDATA\n"
            file.write(header)
            OG2color = dict()
            for index, row in df.iterrows():
                acc, OG = index.split('.')[0], row[rank]
                if OG not in OG2color:
                    color = '#' + "%06x" % random.randint(0, 0xFFFFFF)
                    OG2color[OG] = color
                file.write(f"{acc}\t{OG2color[OG]}\t{OG}\n")
# make_OG_files()
                
# # Make seed label files
# seed_df = pd.read_csv('data/seeds-enriched.tsv', sep='\t')
# make_taxonomy_label_files(seed_df)
# make_activity_label_file(seed_df)
# make_binary_label_files(seed_df)
# make_aguilera_subclass_label_file(seed_df)
# # make_aguilera_subclass_label_file_text(seed_df)
# make_domain_label_file(seed_df)

# # Make blast hit label files
# df_blast_hits = pd.read_csv('data/blast/unique-hits-enriched-interproscan.tsv', sep='\t')
# make_domain_label_file(df_blast_hits, blast_hits=True)
# make_taxonomy_label_files(df_blast_hits, blast_hits=True)

# # Make aguilera label files
# df_aguilera = pd.read_excel('data/Aguilera-data/aguilera_with_seq.xlsx')
# make_aguilera_subclass_label_file_text(df_aguilera)

# Uniprot
df_uniprot_hits = pd.read_csv('data/pfam/protein-matching-PF00264-interproscan2.tsv', sep='\t')
# make_domain_label_file(df_uniprot_hits, uniprot_hits=True)
make_taxonomy_label_files(df_uniprot_hits)
# make_score_label_file()
# make_coverage_label_file()
# make_match_length_file()
make_taxonomy_arrow_files(df_uniprot_hits, 'all')

# make_taxonomy_arrow_files(df_uniprot_hits, 'order', 'all')
# make_taxonomy_arrow_files(df_uniprot_hits, 'class', 'all')
# make_taxonomy_arrow_files(df_uniprot_hits, 'order', 'fungal')
# make_taxonomy_arrow_files(df_uniprot_hits, 'family', 'fungal')
# make_number_of_copies_file()

# Species tree
# df_species_tree = pd.read_csv('species.tsv', sep='\t')
# make_taxonomy_files_species_tree(df_species_tree)