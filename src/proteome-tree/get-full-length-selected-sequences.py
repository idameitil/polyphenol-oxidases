from Bio import SeqIO

def get_fasta_sequences(fasta_filename):
    fasta_alignment = SeqIO.parse(fasta_filename, 'fasta')
    acc2seq = {}
    for fasta in fasta_alignment:
        acc2seq[fasta.id.split('|')[0]] = fasta.seq
    return acc2seq

acc2seq_full_length = get_fasta_sequences('data/pfam/protein-matching-PF00264.fasta')

def get_selected_ids(domain, rank):
    fasta_filename = f'data/proteome-tree/{domain}-one_proteome_per_{rank}.fa'
    fasta_selected = SeqIO.parse(fasta_filename, 'fasta')
    ids_selected = [entry.id for entry in fasta_selected]
    return ids_selected

def write_full_length_fasta(domain, rank):
    ids_selected = get_selected_ids(domain, rank)
    output_filename = f"data/proteome-tree/{domain}-one_proteome_per_{rank}.full-length.fa"
    with open(output_filename, 'w') as outfile:
        for acc in acc2seq_full_length:
            if acc in ids_selected:
                outfile.write(f">{acc}\n{acc2seq_full_length[acc]}\n")

write_full_length_fasta('all', 'order')
