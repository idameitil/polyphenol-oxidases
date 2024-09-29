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


entries = [
    {
        "pdb_name": "2y9x",
        "alignment_name": "TtPPO",
        "family": "a",
        "family_name": "a",
        "printed_name": "i1, TtPPO",
        "descriptive_name": "VvPPO",
        "domain_start_structure": 77,
        "domain_end": 308,
        "threshold": 0.99,
    },
    {
        "pdb_name": "2y9x",
        "alignment_name": "AoCO4",
        "family": "b",
        "family_name": "b",
        "printed_name": "i2, AoCO4",
        "descriptive_name": "EdHc",
        "domain_start_structure": 52,
        "domain_end": 308,
        "threshold": 0.99,
    },
    {
        "pdb_name": "2y9x",
        "alignment_name": "CgPPO-266",
        "family": "c",
        "family_name": "c",
        "printed_name": "i3, CgPPO-266",
        "descriptive_name": "AbTyr",
        "domain_start_structure": 52,
        "domain_end": 308,
        "threshold": 0.99,
    },
    {
        "pdb_name": "2y9x",
        "alignment_name": "MtPPO7",
        "family": "d",
        "family_name": "d",
        "printed_name": "i4, MtPPO7",
        "descriptive_name": "A0A7M6DMJ9",
        "domain_start_structure": 52,
        "domain_end": 308,
        "threshold": 0.99,
    },
    {
        "pdb_name": "2y9x",
        "alignment_name": "MtPPO-809",
        "family": "e",
        "family_name": "e",
        "printed_name": "i5, MtPPO-809",
        "descriptive_name": "A0A7M6DMJ9",
        "domain_start_structure": 52,
        "domain_end": 308,
        "threshold": 0.99,
    },
]

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
alignment_all_groups = "data/compare-architectures/aligned-sequences-noallgaps.fa"
alignment_all_dict = read_MSA_file(alignment_all_groups)

# Get conserved
for i in range(len(entries)):
    entry = entries[i]
    clade_dir = "data/mrbayes/short-fungal-0507/clades"
    alignment_filename = f"{clade_dir}/alignments/{entry['family']}-linsi.fa"
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

# Write conserved file
outfilename = f"data/compare-architectures/fungi/conserved_{threshold}.js"
with open(outfilename, "w") as outfile:
    outfile.write("const conservedResidues = {\n")
    for i in range(len(entries)):
        entry = entries[i]
        outfile.write(f"\t'{entry['printed_name']}': {entry['conserved_positions']},\n")
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
outfilename = f"data/compare-architectures/fungi/architecture.js"
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
outfilename = f"data/compare-architectures/fungi/lengths.js"
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
