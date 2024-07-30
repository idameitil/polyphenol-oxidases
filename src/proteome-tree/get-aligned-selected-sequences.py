from Bio import SeqIO
import pandas as pd
import sys

def get_fasta_sequences(fasta_filename):
    fasta = SeqIO.parse(fasta_filename, "fasta")
    acc2seq = {}
    for fasta in fasta:
        acc = fasta.id.split("|")[0]
        acc2seq[acc] = fasta.seq
    return acc2seq


def get_fasta_sequences_seeds(fasta_filename):
    fasta_alignment = SeqIO.parse(fasta_filename, "fasta")
    acc2seq = {}
    for fasta in fasta_alignment:
        acc2seq[fasta.id] = fasta.seq
    return acc2seq


acc2seq = get_fasta_sequences("data/pfam/protein-matching-PF00264.fasta")
acc2seq_seeds_trimmed = get_fasta_sequences_seeds("data/seeds-trimmed.fa")
seeds_include = [
    "ScTyr",
    "BmTyr",
    "VvPPO",
    "IbCO",
    "EdHc",
    "RvHc",
    "CgPPO-473",
    "CgPPO-266",
    "MtPPO-809",
    "MtPPO-010",
    "PpPPO-c2092",
    "MtPPO7",
    "AbTyr",
    "NpsF",
    "SlPPO1",
    "AoMelO",
    "HjTyr",
    "PsTyr",
]

def read_ids(filename):
    id2position = {}
    with open(filename) as infile:
        selected = [line.strip() for line in infile]
    return selected

def write_fasta(outfilename, included, seeds_include=[]):
    with open(outfilename, 'w') as outfile:
        for entry in included:
            acc = entry.split(',')[0]
            start = int(entry.split(',')[1].split('-')[0])
            stop = int(entry.split(',')[1].split('-')[1])
            trimmed_seq = acc2seq[acc][start:stop]
            outfile.write(f">{acc}/{start}-{stop}\n{trimmed_seq}\n")
        for acc in acc2seq_seeds_trimmed:
            if acc in seeds_include:
                outfile.write(f">{acc}\n{acc2seq_seeds_trimmed[acc]}\n")

def write_noOverlap_filtered_fasta():
    selected = read_ids("data/proteome-tree/sequenceIds-fungi-order-filtered-noOverlap.txt")
    output_filename = f"data/proteome-tree/sequences-fungi-order-filtered-noOverlap.trimmed.fa"
    write_fasta(output_filename, included=selected)

def write_noOverlap_removed_fasta():
    selected = read_ids("data/proteome-tree/sequenceIds-fungi-order-removed-noOverlap.txt")
    output_filename = f"data/proteome-tree/sequences-fungi-order-removed-noOverlap.trimmed.fa"
    write_fasta(output_filename, included=selected)

def write_filtered_sequences_fasta(domain, rank):
    selected = read_ids(f"data/proteome-tree/sequenceIds-{domain}-{rank}-filtered.txt")
    output_filename = f"data/proteome-tree/sequences-{domain}-{rank}-filtered.trimmed.fa"
    if domain == 'all':
        write_fasta(output_filename, included=selected, seeds_include=seeds_include)
    else:
        write_fasta(output_filename, included=selected)

def write_removed_sequences_fasta(domain, rank):
    filtered_out_ids = read_ids(f"data/proteome-tree/sequenceIds-{domain}-{rank}-removed.txt")
    output_filename = f"data/proteome-tree/sequences-{domain}-{rank}-removed.trimmed.fa"
    write_fasta(outfilename=output_filename, included=filtered_out_ids)

def write_all_sequences_fasta(domain, rank):
    selected = read_ids(f"data/proteome-tree/sequenceIds-{domain}-{rank}-filtered.txt")
    filtered_out_ids = read_ids(f"data/proteome-tree/sequenceIds-{domain}-{rank}-removed.txt")
    included = selected + filtered_out_ids
    output_filename = f"data/proteome-tree/sequences-{domain}-{rank}-notFiltered.trimmed.fa"
    write_fasta(outfilename=output_filename, included=included)

write_filtered_sequences_fasta("all", "class")
write_removed_sequences_fasta("all", "class")

write_filtered_sequences_fasta("fungi", "order")
write_removed_sequences_fasta("fungi", "order")
write_noOverlap_filtered_fasta()
write_noOverlap_removed_fasta()
write_all_sequences_fasta('fungi', 'order')