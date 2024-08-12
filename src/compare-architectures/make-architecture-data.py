import sys, os
import pandas as pd
import json
from Bio import SeqIO
import numpy as np
from scipy import stats

threshold = 0.95
aminoacids = [
    "A",
    "C",
    "D",
    "E",
    "F",
    "G",
    "H",
    "I",
    "K",
    "L",
    "M",
    "N",
    "P",
    "Q",
    "R",
    "S",
    "T",
    "V",
    "W",
    "Y",
    "-",
    "X",
    "J",
]


def get_name(entry):
    return f"{entry['family_name']}, {entry['descriptive_name']}"


def read_MSA_file(MSA_filename):
    with open(MSA_filename, "r") as MSA_file:
        proteins = SeqIO.parse(MSA_file, "fasta")
        fasta_dict = {protein.id: protein.seq for protein in proteins}
    return fasta_dict


def AA_to_number(AA):
    try:
        return aminoacids.index(AA)
    except:
        return aminoacids.index("X")


def number_to_AA(number):
    return aminoacids[number]


AA_to_number_vectorized = np.vectorize(AA_to_number)
number_to_AA_vectorized = np.vectorize(number_to_AA)


def get_conserved_residues(fasta_dict, threshold=0.95, include_aliphatic=True):
    if include_aliphatic:
        AAs_ignore = ["-"]
    else:
        # AAs_ignore = ['-', 'G', 'A', 'V', 'C', 'P', 'L', 'I', 'M', 'W', 'F']
        AAs_ignore = ["-", "G", "A", "V", "L", "I", "M"]
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
    return {
        pos: {"AA": AA, "freq": freq}
        for pos, AA, freq in zip(
            conserved_positions, conserved_AAs, frequencies_conserved
        )
    }


def get_specific_positions_conserved_residues(
    accession, conserved_residues, fasta_dict
):
    """Gets positions of conserved residues in a specific protein"""
    seq = fasta_dict[accession]
    protein_position = 0
    positions = []
    alignment_position = 0
    for AA in seq:
        if AA != "-":
            protein_position += 1
        if alignment_position in conserved_residues:
            if AA == conserved_residues[alignment_position]["AA"]:
                positions.append(
                    {
                        "pos": protein_position,
                        "AA": AA,
                        "freq": conserved_residues[alignment_position]["freq"],
                    }
                )
            else:
                print(
                    f"Warning: {accession} doesn't have the conserved residue {alignment_position} {conserved_residues[alignment_position]}"
                )
        alignment_position += 1
    return positions


# entries = [
#     {'acc': '2p3x', 'family': 'a_plants', 'family_name': 'a (plant)', 'descriptive_name': 'VvPPO', 'domain_start_structure': 77, 'domain_end': 284, 'threshold':0.9},
#     # {'acc': 'Q63JI5', 'family': 'b_cnidaria', 'descriptive_name': 'Q63JI5', 'domain_start_structure': 49, 'domain_end': 334, 'threshold':0.9},
#     {'acc': '2y9x', 'family': 'c_long_fungal', 'family_name': 'c (long fungal)', 'descriptive_name': 'AbPPO3', 'domain_start_structure': 52, 'domain_end': 308, 'threshold':0.95},
#     {'acc': '1wx2', 'family': 'd_bacteria', 'family_name': 'd (bacteria)', 'descriptive_name': 'ScTyr', 'domain_start_structure': 29, 'domain_end': 229, 'threshold':0.95},
#     # {'acc': '1wx2', 'family': 'k_bacteria2', 'descriptive_name': 'ScTyr', 'domain_start_structure': 29, 'domain_end': 226, 'threshold':0.95},
#     {'acc': '5m8l', 'family': 'e_chordata', 'family_name': 'e (chordata)', 'descriptive_name': 'TyrHs', 'domain_start_structure': 184, 'domain_end': 416, 'threshold':0.95},
#     {'acc': 'V3ZAB2', 'family': 'f_mollusc', 'family_name': 'f (mollusc)', 'descriptive_name': 'V3ZAB2', 'domain_start_structure': 130, 'domain_end': 305, 'threshold':0.93},
#     {'acc': 'D0N318', 'family': 'h_oomycota', 'family_name': 'h (oomycota)', 'descriptive_name': 'D0N318', 'domain_start_structure': 77, 'domain_end': 279, 'threshold':0.99},
#     {'acc': '4j3p', 'family': 'i_short_fungal', 'family_name': 'i (short fungal)', 'descriptive_name': 'AoCO4', 'domain_start_structure': 93, 'domain_end': 324, 'threshold':0.88},
#     {'acc': 'F4PFF7', 'family': 'j_zoopagomycota', 'family_name': 'j (zoopagomycota)', 'descriptive_name': 'F4PFF7', 'domain_start_structure': 68, 'domain_end': 251, 'threshold':0.92},
# ]
entries = [
    {
        "pdb_name": "2p3x",
        "alignment_name": "VvPPO",
        "family": "a_plants",
        "family_name": "a",
        "printed_name": "a, 2P3X (VvPPO)",
        "descriptive_name": "VvPPO",
        "domain_start_structure": 77,
        "domain_end": 284,
        "threshold": 0.95,
    },
    {
        "pdb_name": "1lnl",
        "alignment_name": "RvHc",
        "family": "b_hc",
        "family_name": "b",
        "printed_name": "b, 1LNL (RvHc)",
        "descriptive_name": "RvHc",
        "domain_start_structure": 32,
        "domain_end": 224,
        "threshold": 0.99,
    },
    {
        "pdb_name": "2y9x",
        "alignment_name": "AbTyr",
        "family": "c_long_fungal",
        "family_name": "c",
        "printed_name": "c, 2Y9X (AbTyr)",
        "descriptive_name": "AbTyr",
        "domain_start_structure": 52,
        "domain_end": 308,
        "threshold": 0.95,
    },
    {
        "pdb_name": "5m8l",
        "alignment_name": "P17643.20184-416",
        "family": "d_chordata",
        "family_name": "d",
        "printed_name": "d, 5M8L (HsTrp1)",
        "descriptive_name": "HsTrp1",
        "domain_start_structure": 184,
        "domain_end": 416,
        "threshold": 0.95,
    },
    {
        "pdb_name": "V3ZAB2",
        "alignment_name": "V3ZAB2.10130-305",
        "family": "e_mollusc",
        "family_name": "e",
        "printed_name": "e, V3ZAB2",
        "descriptive_name": "V3ZAB2",
        "domain_start_structure": 130,
        "domain_end": 305,
        "threshold": 0.95,
    },
    {
        "pdb_name": "A0A7M6DMJ9",
        "alignment_name": "A0A7M6DMJ9.10298-479",
        "family": "f_cnidaria",
        "family_name": "f",
        "printed_name": "f, A0A7M6DMJ9",
        "descriptive_name": "A0A7M6DMJ9",
        "domain_start_structure": 298,
        "domain_end": 479,
        "threshold": 0.95,
    },
    {
        "pdb_name": "D0N318",
        "alignment_name": "D0N318.1077-279",
        "family": "g_oomycota",
        "family_name": "g",
        "printed_name": "g, D0N318",
        "descriptive_name": "D0N318",
        "domain_start_structure": 77,
        "domain_end": 279,
        "threshold": 0.95,
    },
    {
        "pdb_name": "4j3p",
        "alignment_name": "A0A1S9DJX3.10118-349",
        "family": "h_short_fungal",
        "family_name": "h",
        "printed_name": "h, 4J3P (AoCO4)",
        "descriptive_name": "AoCO4",
        "domain_start_structure": 93,
        "domain_end": 324,
        # "threshold": 0.88,
        "threshold": 0.90,
    },
    {
        "pdb_name": "F4PFF7",
        "alignment_name": "F4PFF7.1068-251",
        "family": "i_zoopago",
        "family_name": "i",
        "printed_name": "i, F4PFF7",
        "descriptive_name": "F4PFF7",
        "domain_start_structure": 68,
        "domain_end": 251,
        "threshold": 0.95,
    },
]


# A0A1S9DJX3.10118-349
def make_sequencePos2alignedPos(aligned):
    not_aligned_counter = 0
    aligned_counter = 0
    sequencePos2alignedPos = {}
    for char in aligned:
        aligned_counter += 1
        if char != "-":
            not_aligned_counter += 1
            sequencePos2alignedPos[not_aligned_counter] = aligned_counter
    return sequencePos2alignedPos


# Read alignment all
# alignment_all_groups = 'data/compare-architectures/aligned-sequences-noallgaps.fa'
# alignment_all_groups = 'data/compare-architectures/linsi2-noallgaps.fa'
alignment_all_groups = "data/compare-architectures/aligned-sequences-noallgaps.fa"
alignment_all_dict = read_MSA_file(alignment_all_groups)

# Get conserved
for i in range(len(entries)):
    entry = entries[i]
    # clade_dir = "data/mrbayes/all/clades"
    clade_dir = "data/mrbayes/all-seeds-0619/clades"
    # alignment_filename = f"{clade_dir}/alignments/{entry['family']}-linsi.fa"
    alignment_filename = f"{clade_dir}/alignments-nofragments/{entry['family']}.fa"
    fasta_dict = read_MSA_file(alignment_filename)
    conserved_residues = get_conserved_residues(
        fasta_dict, threshold=entry["threshold"], include_aliphatic=True
    )
    positions = get_specific_positions_conserved_residues(
        entry["alignment_name"], conserved_residues, fasta_dict
    )
    entries[i]["conserved_positions"] = {
        position["pos"] + entry["domain_start_structure"] - 1: position["AA"]
        for position in positions
    }
    # Get positions in alignment
    sequencePos2alignedPos = make_sequencePos2alignedPos(
        alignment_all_dict[entries[i]["alignment_name"]]
    )
    entries[i]["conserved_positions_alignment"] = {
        sequencePos2alignedPos[position["pos"]]
        + entry["domain_start_structure"]
        - 1: position["AA"]
        for position in positions
    }

# Write conserved file
outfilename = f"data/compare-architectures/conserved_{threshold}.js"
with open(outfilename, "w") as outfile:
    outfile.write("const conservedResidues = {\n")
    for i in range(len(entries)):
        entry = entries[i]
        # outfile.write(f"\t'{get_name(entry)}': {entry['conserved_positions']},\n")
        outfile.write(f"\t'{entry['printed_name']}': {entry['conserved_positions']},\n")
    outfile.write("};")

# Write conserved file - alignment positions
outfilename = f"data/compare-architectures/conserved_{threshold}-alignment-positions.js"
with open(outfilename, "w") as outfile:
    outfile.write("const conservedResidues = {\n")
    for i in range(len(entries)):
        entry = entries[i]
        outfile.write(
            f"\t'{entry['printed_name']}': {entry['conserved_positions_alignment']},\n"
        )
    outfile.write("};")

# Get architecture strings
table_folder = "data/compare-architectures/architecture-tables-dssp"
for i in range(len(entries)):
    entry = entries[i]
    table_filename = f"{table_folder}/{entry['pdb_name']}.csv"
    df = pd.read_csv(table_filename, sep=",")
    architecture_string = ""
    previous = 0
    for index, row in df.iterrows():
        length = int(row.end) - previous
        architecture_string += row.type * length
        previous = int(row.end)
    entries[i]["architecture_string"] = architecture_string[
        entry["domain_start_structure"] - 1 : entry["domain_end"]
    ]
    entries[i]["length"] = entry["domain_end"] - entry["domain_start_structure"]

# Write architecture file
outfilename = f"data/compare-architectures/architecture.js"
with open(outfilename, "w") as outfile:
    outfile.write("const architectures = {\n")
    for i in range(len(entries)):
        entry = entries[i]
        outfile.write(
            f"\t'{entry['printed_name']}': {{\n\t\t'domain_start_structure': {entry['domain_start_structure']},\n\t\t'architecture_string': '{entry['architecture_string']}'\n\t\t}}"
        )
        if i < len(entries) - 1:
            outfile.write(",\n")
        else:
            outfile.write("\n")
    outfile.write("};")

# Write length file
outfilename = f"data/compare-architectures/lengths.js"
with open(outfilename, "w") as outfile:
    outfile.write("const lengths = {\n")
    max_length = 0
    for i in range(len(entries)):
        entry = entries[i]
        outfile.write(f"\t'{entry['printed_name']}': {entry['length']},\n")
        if entry["length"] > max_length:
            max_length = entry["length"]
    outfile.write("};\n")
    outfile.write(f"const max_length = {max_length};\n")
