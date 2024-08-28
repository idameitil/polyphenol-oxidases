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

# Uniprot accessions of seeds that are in the genome dataset
uniprot_seeds = {
    "P17643": "HsTrp1",
    "A0A1S9DJX3": "AoCO4",
    "C7FF04": "AbPPO3",
    "C7FF05": "AbPPO4",
    "A0A8H7C4Z8": "AbPPO6",
    "A0A8H7C551": "AbPPO1",
    "A0A8H7F0X6": "AbPPO3",
    "A0A8H7F114": "AbPPO2",
    "A0A8H7F118": "AbPPO4",
    "A0A8H7F178": "AbPPO5",
    "A0A1S9DK56": "AoMelB",
    "A0A498IXV7": "MdPPO1",
    "D5DZK6": "BmTyr"
}

fungal_seeds = [
    'AoMelB',
    'AoCO4',
    'CgPPO-473',
    'CgPPO-266',
    'MtPPO-809',
    'MtPPO-010',
    'PpPPO-c2092',
    'MtPPO7',
    'AbPPO1',
    'AbPPO2',
    'AbPPO3',
    'AbPPO4',
    'AbPPO5',
    'AbPPO6',
    'HjTyr',
    'PsTyr'
]

seeds_exclude = ['A0A6P7SER3']

def read_ids(filename):
    with open(filename) as infile:
        selected = [line.strip() for line in infile]
    return selected


def write_fasta(outfilename, included, include_seeds=False, include_fungal_seeds=False):
    with open(outfilename, "w") as outfile:
        for entry in included:
            acc = entry.split(",")[0]
            if include_seeds or include_fungal_seeds:
                if acc in uniprot_seeds:
                    continue
            if acc in seeds_exclude:
                continue
            start = int(entry.split(",")[1].split("-")[0])
            stop = int(entry.split(",")[1].split("-")[1])
            trimmed_seq = acc2seq[acc][start:stop]
            outfile.write(f">{acc}/{start}-{stop}\n{trimmed_seq}\n")
        if include_seeds:
            for acc in acc2seq_seeds_trimmed:
                outfile.write(f">{acc}\n{acc2seq_seeds_trimmed[acc]}\n")
        if include_fungal_seeds:
            for acc in acc2seq_seeds_trimmed:
                if acc in fungal_seeds:
                    outfile.write(f">{acc}\n{acc2seq_seeds_trimmed[acc]}\n") 
            

def write_noOverlap_filtered_fasta(lignin_degraders = False):
    if not lignin_degraders:
        selected = read_ids("data/proteome-tree/sequenceIds-fungi-order-filtered-noOverlap.txt")
        output_filename = (f"data/proteome-tree/sequences-fungi-order-filtered-noOverlap.trimmed.fa")
    else:
        selected = read_ids("data/proteome-tree/sequenceIds-fungi-with-lignin-degraders-order-filtered-noOverlap.txt")
        output_filename = (f"data/proteome-tree/sequences-fungi-with-lignin-degraders-order-filtered-noOverlap.trimmed.fa")
    write_fasta(output_filename, included=selected)


def write_noOverlap_removed_fasta(lignin_degraders = False):
    if not lignin_degraders:
        selected = read_ids("data/proteome-tree/sequenceIds-fungi-order-removed-noOverlap.txt")
        output_filename = (f"data/proteome-tree/sequences-fungi-order-removed-noOverlap.trimmed.fa")
    else:
        selected = read_ids("data/proteome-tree/sequenceIds-fungi-with-lignin-degraders-order-removed-noOverlap.txt")
        output_filename = (f"data/proteome-tree/sequences-fungi-with-lignin-degraders-order-removed-noOverlap.trimmed.fa")
    write_fasta(output_filename, included=selected)


def write_filtered_sequences_fasta(domain, rank):
    selected = read_ids(f"data/proteome-tree/sequenceIds-{domain}-{rank}-filtered.txt")
    output_filename = (f"data/proteome-tree/sequences-{domain}-{rank}-filtered.trimmed.fa")
    write_fasta(output_filename, included=selected)
    if domain == "all":
        output_filename = (f"data/proteome-tree/sequences-{domain}-{rank}-filtered-andseeds.trimmed.fa")
        write_fasta(output_filename, included=selected, include_seeds=True)
    else:
        output_filename = (f"data/proteome-tree/sequences-{domain}-{rank}-filtered-andseeds.trimmed.fa")
        write_fasta(output_filename, included=selected, include_fungal_seeds=True)

def write_removed_sequences_fasta(domain, rank):
    filtered_out_ids = read_ids(f"data/proteome-tree/sequenceIds-{domain}-{rank}-removed.txt")
    output_filename = f"data/proteome-tree/sequences-{domain}-{rank}-removed.trimmed.fa"
    write_fasta(outfilename=output_filename, included=filtered_out_ids)


def write_all_sequences_fasta(domain, rank):
    # selected = read_ids(f"data/proteome-tree/sequenceIds-{domain}-{rank}-filtered.txt")
    # filtered_out_ids = read_ids(f"data/proteome-tree/sequenceIds-{domain}-{rank}-removed.txt")
    # included = selected + filtered_out_ids
    included = read_ids(f"data/proteome-tree/sequenceIds-{domain}-{rank}-notFiltered.txt")
    output_filename = (f"data/proteome-tree/sequences-{domain}-{rank}-notFiltered.trimmed.fa")
    write_fasta(outfilename=output_filename, included=included)


write_filtered_sequences_fasta("all", "class")
write_removed_sequences_fasta("all", "class")

write_filtered_sequences_fasta("fungi", "order")
write_removed_sequences_fasta("fungi", "order")
write_noOverlap_filtered_fasta()
write_noOverlap_removed_fasta()
write_all_sequences_fasta("fungi", "order")

# write_all_sequences_fasta("fungi-with-lignin-degraders", "order")

write_noOverlap_filtered_fasta(
    lignin_degraders=True
)
write_noOverlap_removed_fasta(
    lignin_degraders=True
)