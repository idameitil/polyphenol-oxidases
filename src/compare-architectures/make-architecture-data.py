import sys, os
import pandas as pd
import json
from Bio import SeqIO
import numpy as np
from scipy import stats

threshold = 0.95

def get_name(entry):
    return f"{entry['acc']} {entry['family']}"

def read_MSA_file(MSA_filename):
    with open(MSA_filename, 'r') as MSA_file:
        proteins = SeqIO.parse(MSA_file, 'fasta')
        fasta_dict = {protein.id: protein.seq for protein in proteins}
    return fasta_dict

aminoacids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
        'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '-', 'X', 'J']

def AA_to_number(AA):
    try:
        return aminoacids.index(AA)
    except:
        return aminoacids.index('X')

def number_to_AA(number):
    return aminoacids[number]

AA_to_number_vectorized = np.vectorize(AA_to_number)
number_to_AA_vectorized = np.vectorize(number_to_AA)

def get_conserved_residues(fasta_dict, threshold=0.95, include_aliphatic=False):
    if include_aliphatic:
        AAs_ignore = ['-']
    else:
        AAs_ignore = ['-', 'G', 'A', 'V', 'C', 'P', 'L', 'I', 'M', 'W', 'F']
    sequences = np.array([np.array(list(fasta_dict[acc])) for acc in fasta_dict])
    no_sequences = sequences.shape[0]
    sequences_numerical = AA_to_number_vectorized(sequences)
    mode = stats.mode(sequences_numerical)
    mode_AAs = number_to_AA_vectorized(mode[0][0])
    frequencies = mode[1][0] / no_sequences
    condition = (frequencies > threshold) & (np.isin(mode_AAs, AAs_ignore, invert=True))
    conserved_AAs = mode_AAs[condition]
    frequencies_conserved = frequencies[condition]
    conserved_positions = list(np.where(condition)[0])
    return {pos: {'AA': AA, 'freq': freq} for pos, AA, freq in zip(conserved_positions, conserved_AAs, frequencies_conserved)}

def get_specific_positions_conserved_residues(accession, conserved_residues, fasta_dict):
    """Gets positions of conserved residues in a specific protein"""
    seq = fasta_dict[accession]
    protein_position = 0
    positions = []
    alignment_position = 0
    for AA in seq:
        if AA != '-':
            protein_position += 1
        if alignment_position in conserved_residues:
            if AA == conserved_residues[alignment_position]['AA']:
                positions.append({'pos':protein_position, 'AA': AA, 'freq': conserved_residues[alignment_position]['freq']})
            else:
                print(f"Warning: {accession} doesn't have the conserved residue {alignment_position} {conserved_residues[alignment_position]}")
        alignment_position += 1
    return positions

entries = [
    {'acc': '4j3p', 'family': 'i_short_fungal', 'descriptive_name': 'AoCO4', 'domain_start': 93, 'domain_end': 322},
    {'acc': '5m8l', 'family': 'e_chordata', 'descriptive_name': 'TyrHs', 'domain_start': 184, 'domain_end': 414},
    {'acc': '1wx2', 'family': 'd_bacteria', 'descriptive_name': 'TyrSc', 'domain_start': 29, 'domain_end': 226}, 
    {'acc': '2y9x', 'family': 'c_long_fungal', 'descriptive_name': 'AbPPO3', 'domain_start': 52, 'domain_end': 308}, 
    {'acc': '2p3x', 'family': 'a_plants', 'descriptive_name': 'PPOVv', 'domain_start': 77, 'domain_end': 284} 
]

# Write manual conservation file
# filename = "data/compare-architectures/conserved_manual.json"
# with open(filename) as infile:
#     conserved_dict = json.load(infile)

# outfilename = f"data/compare-architectures/conserved_manual.js"
# with open(outfilename, 'w') as outfile:
#     outfile.write("const conservedResidues = {\n")
#     for i in range(len(entries)):
#         for j in range(len(entries[i])):
#             entry = entries[i][j]
#             outfile.write(f"\t'{get_name(entry)}': {conserved_dict[entry['acc']]},\n")
#     outfile.write('};')

# Get conserved
for i in range(len(entries)):
    entry = entries[i]
    alignment_filename = f"data/mrbayes/all/clades/alignments/{entry['family']}-linsi.fa"
    fasta_dict = read_MSA_file(alignment_filename)
    conserved_residues = get_conserved_residues(fasta_dict, threshold=threshold, include_aliphatic=True)
    if entry['family'] == 'e_chordata':
        print(conserved_residues)
    positions = get_specific_positions_conserved_residues(entry['descriptive_name'], conserved_residues, fasta_dict)
    # entries[i]['conserved_positions'] = {position['pos'] + entry['domain_start']: position['AA'] for position in positions}
    entries[i]['conserved_positions'] = {position['pos']: position['AA'] for position in positions}

# Write conserved file
outfilename = f"data/compare-architectures/conserved_{threshold}.js"
with open(outfilename, 'w') as outfile:
    outfile.write("const conservedResidues = {\n")
    for i in range(len(entries)):
        entry = entries[i]
        outfile.write(f"\t'{get_name(entry)}': {entry['conserved_positions']},\n")
    outfile.write('};')

# Get architecture strings
table_folder = "data/compare-architectures/architecture-tables-dssp"
for i in range(len(entries)):
    entry = entries[i]
    table_filename = f"{table_folder}/{entry['acc']}.csv"
    df = pd.read_csv(table_filename, sep=',')
    architecture_string = ''
    previous = 0
    for index, row in df.iterrows():
        length = int(row.end) - previous
        architecture_string += row.type * length
        previous = int(row.end)
    entries[i]['architecture_string'] = architecture_string[entry['domain_start']:entry['domain_end']]
    entries[i]['length'] = entry['domain_end'] - entry['domain_start']

# Write architecture file
outfilename = f"data/compare-architectures/architecture.js"
with open(outfilename, 'w') as outfile:
    outfile.write("const architectures = {\n")
    for i in range(len(entries)):
        entry = entries[i]
        outfile.write(f"\t'{get_name(entry)}': '{entry['architecture_string']}'")
        if i < len(entries) - 1:
            outfile.write(",\n")
        else:
            outfile.write("\n")
    outfile.write('};')

# Write length file
outfilename = f"data/compare-architectures/lengths.js"
with open(outfilename, 'w') as outfile:
    outfile.write("const lengths = {\n")
    max_length = 0
    for i in range(len(entries)):
        entry = entries[i]
        outfile.write(f"\t'{get_name(entry)}': {entry['length']},\n")
        if entry['length'] > max_length:
            max_length = entry['length']
    outfile.write("};\n")
    outfile.write(f"const max_length = {max_length};\n")

# accessions_filename = "data/garry-figure/family-representatives.tsv"
# df = pd.read_csv(accessions_filename, sep='\t', names=['family', 'accession'])
# threshold = 0.95

# def replacer(s, newstring, index, nofail=False):
#     # raise an error if index is outside of the string
#     if not nofail and index not in range(len(s)):
#         raise ValueError("index outside given string")

#     # if not erroring, but the index is still not in the correct range..
#     if index < 0:  # add it to the beginning
#         return newstring + s
#     if index > len(s):  # add it to the end
#         return s + newstring

#     # insert the new string between "slices" of the original
#     return s[:index] + newstring + s[index + 1:]

# outfilename = "data/garry-figure/conserved.out"
# with open(outfilename, 'w') as outfile:
#     for index, row in df.iterrows():
#         outfile.write(row.family + '\n')
#         outfile.write(row.accession + '\n')
#         alignment_filename = f"data/hhblits_cazy_families/msas-family-names/{row.family}.fa"
#         fasta_dict = read_MSA_file(alignment_filename)
#         conserved_residues = get_conserved_residues(fasta_dict, threshold=threshold)
#         positions = get_specific_positions_conserved_residues(row.accession, conserved_residues, fasta_dict)

#         for position in positions:
#             strings[row.accession] = replacer(s = strings[row.accession], newstring = position['AA'], index = position['pos'])

# with open(outfilename_retaining, 'w') as outfile:
#     for family in retaining:
#         acc = family2acc[family]
#         outfile.write(f"{family}, {acc}\n")
#         outfile.write(strings[acc] + '\n')