from Bio import SeqIO

def get_fasta_sequences(fasta_filename):
    fasta_alignment = SeqIO.parse(fasta_filename, 'fasta')
    acc2seq = {}
    for fasta in fasta_alignment:
        acc2seq[fasta.id.split('.')[0]] = fasta.seq
    return acc2seq

acc2seqhmmalign = get_fasta_sequences('data/pfam/PF00264.alignment.uniprot-cleaned.fa')
acc2seqtrimmed = get_fasta_sequences('data/pfam/PF00264.alignment.uniprot-nogaps.fa')

def get_selected_ids(domain, rank):
    fasta_filename = f'data/proteome-tree/{domain}-one_proteome_per_{rank}.fa'
    fasta_selected = SeqIO.parse(fasta_filename, 'fasta')
    ids_selected = [entry.id for entry in fasta_selected]
    return ids_selected

def write_hmmalign_fasta(domain, rank):
    ids_selected = get_selected_ids(domain, rank)
    output_filename = f"data/proteome-tree/{domain}-one_proteome_per_{rank}.hmmalign.fa"
    with open(output_filename, 'w') as outfile:
        for acc in acc2seqhmmalign:
            if acc in ids_selected:
                outfile.write(f">{acc}\n{acc2seqhmmalign[acc]}\n")

def write_trimmed_fasta(domain, rank):
    """for mafft tree"""
    ids_selected = get_selected_ids(domain, rank)
    output_filename = f"data/proteome-tree/{domain}-one_proteome_per_{rank}.trimmed.fa"
    with open(output_filename, 'w') as outfile:
        for acc in acc2seqtrimmed:
            if acc in ids_selected:
                outfile.write(f">{acc}\n{acc2seqtrimmed[acc]}\n")

write_hmmalign_fasta('fungal', 'order')
write_hmmalign_fasta('fungal', 'family')
write_hmmalign_fasta('all', 'class')
write_hmmalign_fasta('all', 'order')

write_trimmed_fasta('fungal', 'order')
write_trimmed_fasta('fungal', 'family')
write_trimmed_fasta('all', 'class')
write_trimmed_fasta('all', 'order')