from Bio import SeqIO
import pandas as pd

def get_fasta_sequences(fasta_filename = 'data/pfam/PF00264.alignment.uniprot.fa'):
    fasta_alignment = SeqIO.parse(fasta_filename, 'fasta')
    acc2seq = {}
    for fasta in fasta_alignment:
        acc = fasta.id.split('.')[0]
        version = fasta.id.split('.')[1][0]
        position = fasta.id.split('/')[1]
        trimmed_seq = fasta.seq.replace('-', '')
        acc2seq[f"{acc},{position}"] = {'acc': acc, 
                                        'version': version, 
                                        'position': position, 
                                        'trimmed_seq': trimmed_seq}
    return acc2seq

acc2seq_trimmed = get_fasta_sequences()

def get_selected_ids():
    infile_name = f'data/proteome-tree/selected-sequences-all-class.txt'
    with open(infile_name) as infile:
        selected = [line.strip() for line in infile]
    return selected

def write_trimmed_fasta(domain, rank):
    """for mafft tree"""
    selected = get_selected_ids()
    output_filename = f"data/proteome-tree/{domain}-one_proteome_per_{rank}.trimmed.fa"
    with open(output_filename, 'w') as outfile:
        for entry in acc2seq_trimmed:
            if entry in selected:
                outfile.write(f">{acc2seq_trimmed[entry]['acc']}.{acc2seq_trimmed[entry]['version']}/{acc2seq_trimmed[entry]['position']}\n{acc2seq_trimmed[entry]['trimmed_seq']}\n")

write_trimmed_fasta('all', 'class')