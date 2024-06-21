from Bio import SeqIO
import pandas as pd


def get_fasta_sequences_location(fasta_filename):
    fasta_alignment = SeqIO.parse(fasta_filename, "fasta")
    acc2seq = {}
    for fasta in fasta_alignment:
        acc = fasta.id.split(".")[0]
        version = fasta.id.split(".")[1][0]
        position = fasta.id.split("/")[1]
        trimmed_seq = fasta.seq.replace("-", "")
        acc2seq[f"{acc},{position}"] = {
            "acc": acc,
            "version": version,
            "position": position,
            "trimmed_seq": trimmed_seq,
        }
    return acc2seq


def get_fasta_sequences_seeds(fasta_filename):
    fasta_alignment = SeqIO.parse(fasta_filename, "fasta")
    acc2seq = {}
    for fasta in fasta_alignment:
        acc2seq[fasta.id] = fasta.seq
    return acc2seq


acc2seq_trimmed = get_fasta_sequences_location("data/pfam/PF00264.alignment.uniprot.fa")
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


def get_selected_ids(domain, rank):
    infile_name = f"data/proteome-tree/selected-sequences-{domain}-{rank}.txt"
    with open(infile_name) as infile:
        selected = [line.strip() for line in infile]
    return selected


def get_filtered_out_ids():
    infile_name = f"data/proteome-tree/filtered-out-sequences-all-class.txt"
    with open(infile_name) as infile:
        filtered_out = [line.strip() for line in infile]
    return filtered_out


def write_trimmed_fasta(domain, rank):
    selected = get_selected_ids(domain, rank)
    output_filename = f"data/proteome-tree/{domain}-one_proteome_per_{rank}.trimmed.fa"
    with open(output_filename, "w") as outfile:
        if domain == 'all':
            for entry in acc2seq_trimmed:
                if entry in selected:
                    outfile.write(
                        f">{acc2seq_trimmed[entry]['acc']}.{acc2seq_trimmed[entry]['version']}/{acc2seq_trimmed[entry]['position']}\n{acc2seq_trimmed[entry]['trimmed_seq']}\n"
                    )
            for acc in acc2seq_seeds_trimmed:
                if acc in seeds_include:
                    outfile.write(f">{acc}\n{acc2seq_seeds_trimmed[acc]}\n")
        if domain == 'fungi':
            all_selected = get_selected_ids(domain='all', rank='class')
            for entry in acc2seq_trimmed:
                # check if it's in the selection of all
                id_formatted = f"{acc2seq_trimmed[entry]['acc']},{acc2seq_trimmed[entry]['position']}"
                if id_formatted in selected and id_formatted not in all_selected:
                    outfile.write(
                        f">{acc2seq_trimmed[entry]['acc']}.{acc2seq_trimmed[entry]['version']}/{acc2seq_trimmed[entry]['position']}\n{acc2seq_trimmed[entry]['trimmed_seq']}\n"
                    ) 


def write_filtered_out():
    filtered_out_ids = get_filtered_out_ids()
    output_filename = f"data/proteome-tree/filtered-out-all-one_proteome_per_class.trimmed.fa"
    with open(output_filename, "w") as outfile:
        for entry in acc2seq_trimmed:
            if entry in filtered_out_ids:
                outfile.write(
                    f">{acc2seq_trimmed[entry]['acc']}.{acc2seq_trimmed[entry]['version']}/{acc2seq_trimmed[entry]['position']}\n{acc2seq_trimmed[entry]['trimmed_seq']}\n"
                )

write_trimmed_fasta("all", "class")
write_trimmed_fasta("fungi", "order")
write_filtered_out()