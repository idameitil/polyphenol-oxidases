import sys, os
import pandas as pd
import json
from Bio import SeqIO
import numpy as np
from scipy import stats

threshold = 0.95
aminoacids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
        'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '-', 'X', 'J']

def get_name(entry):
    return f"{entry['acc']} {entry['family']}"

def read_MSA_file(MSA_filename):
    with open(MSA_filename, 'r') as MSA_file:
        proteins = SeqIO.parse(MSA_file, 'fasta')
        fasta_dict = {protein.id: protein.seq for protein in proteins}
    return fasta_dict

def AA_to_number(AA):
    try:
        return aminoacids.index(AA)
    except:
        return aminoacids.index('X')

def number_to_AA(number):
    return aminoacids[number]

AA_to_number_vectorized = np.vectorize(AA_to_number)
number_to_AA_vectorized = np.vectorize(number_to_AA)

def get_conserved_residues(fasta_dict, threshold=0.95, include_aliphatic=True):
    if include_aliphatic:
        AAs_ignore = ['-']
    else:
        # AAs_ignore = ['-', 'G', 'A', 'V', 'C', 'P', 'L', 'I', 'M', 'W', 'F']
        AAs_ignore = ['-', 'G', 'A', 'V', 'L', 'I', 'M']
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
    {'acc': '2p3x', 'family': 'a_plants', 'descriptive_name': 'PPOVv', 'domain_start_structure': 77, 'domain_end': 284, 'threshold':0.9}, 
    # {'acc': 'Q63JI5', 'family': 'b_cnidaria', 'descriptive_name': 'Q63JI5', 'domain_start_structure': 49, 'domain_end': 334, 'threshold':0.9}, 
    {'acc': '2y9x', 'family': 'c_long_fungal', 'descriptive_name': 'AbPPO3', 'domain_start_structure': 52, 'domain_end': 308, 'threshold':0.95}, 
    {'acc': '1wx2', 'family': 'd_bacteria', 'descriptive_name': 'TyrSc', 'domain_start_structure': 29, 'domain_end': 226, 'threshold':0.95}, 
    # {'acc': '1wx2', 'family': 'k_bacteria2', 'descriptive_name': 'TyrSc', 'domain_start_structure': 29, 'domain_end': 226, 'threshold':0.95}, 
    {'acc': '5m8l', 'family': 'e_chordata', 'descriptive_name': 'TyrHs', 'domain_start_structure': 184, 'domain_end': 416, 'threshold':0.95},
    {'acc': 'V3ZAB2', 'family': 'f_mollusc', 'descriptive_name': 'V3ZAB2', 'domain_start_structure': 130, 'domain_end': 300, 'threshold':0.93},
    {'acc': 'D0N318', 'family': 'h_oomycota', 'descriptive_name': 'D0N318', 'domain_start_structure': 77, 'domain_end': 274, 'threshold':0.99},
    {'acc': '4j3p', 'family': 'i_short_fungal', 'descriptive_name': 'AoCO4', 'domain_start_structure': 93, 'domain_end': 322, 'threshold':0.88},
    {'acc': 'F4PFF7', 'family': 'j_zoopagomycota', 'descriptive_name': 'F4PFF7', 'domain_start_structure': 68, 'domain_end': 247, 'threshold':0.92},
]

# Get conserved
for i in range(len(entries)):
    entry = entries[i]
    alignment_filename = f"data/mrbayes/all/clades/alignments/{entry['family']}-linsi.fa"
    fasta_dict = read_MSA_file(alignment_filename)
    conserved_residues = get_conserved_residues(fasta_dict, threshold=entry['threshold'], include_aliphatic=True)
    positions = get_specific_positions_conserved_residues(entry['descriptive_name'], conserved_residues, fasta_dict)
    entries[i]['conserved_positions'] = {position['pos'] + entry['domain_start_structure'] -1: position['AA'] for position in positions}

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
    entries[i]['architecture_string'] = architecture_string[entry['domain_start_structure']-1:entry['domain_end']-1]
    entries[i]['length'] = entry['domain_end'] - entry['domain_start_structure']

# Write architecture file
outfilename = f"data/compare-architectures/architecture.js"
with open(outfilename, 'w') as outfile:
    outfile.write("const architectures = {\n")
    for i in range(len(entries)):
        entry = entries[i]
        outfile.write(f"\t'{get_name(entry)}': {{\n\t\t'domain_start_structure\': {entry['domain_start_structure']},\n\t\t\'architecture_string\': '{entry['architecture_string']}'\n\t\t}}")
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