from Bio import SeqIO
import pandas as pd

def get_fasta_sequences(fasta_filename):
    fasta_alignment = SeqIO.parse(fasta_filename, 'fasta')
    acc2seq = {}
    for fasta in fasta_alignment:
        acc2seq[fasta.id.split('.')[0]] = fasta.seq
    return acc2seq

acc2seq_hmmalign = get_fasta_sequences('data/pfam/PF00264.alignment.uniprot-cleaned.fa')
acc2seq_trimmed = get_fasta_sequences('data/pfam/PF00264.alignment.uniprot-nogaps.fa')
acc2seq_seeds_hmmalign = get_fasta_sequences('data/seeds-names.hmmalign.fa')
acc2seq_seeds_trimmed = get_fasta_sequences('data/seeds-names.hmmalign-withoutgaps.fa')

def get_selected_ids(domain, rank):
    fasta_filename = f'data/proteome-tree/{domain}-one_proteome_per_{rank}.fa'
    fasta_selected = SeqIO.parse(fasta_filename, 'fasta')
    ids_selected = [entry.id for entry in fasta_selected]
    return ids_selected

def get_fungal_seeds():
    df = pd.read_csv('data/seeds-enriched.tsv', sep='\t')
    return df[df.kingdom == 'Fungi'].descriptive_name.to_list()

fungal_seeds = get_fungal_seeds()

def write_hmmalign_fasta(domain, rank):
    ids_selected = get_selected_ids(domain, rank)
    output_filename = f"data/proteome-tree/{domain}-one_proteome_per_{rank}.hmmalign.fa"
    with open(output_filename, 'w') as outfile:
        for acc in acc2seq_hmmalign:
            if acc in ids_selected:
                outfile.write(f">{acc}\n{acc2seq_hmmalign[acc]}\n")
        for acc in acc2seq_seeds_hmmalign:
            if domain == 'fungal':
                if acc not in fungal_seeds:
                    continue
            outfile.write(f">{acc}\n{acc2seq_seeds_hmmalign[acc]}\n")

def write_trimmed_fasta(domain, rank):
    """for mafft tree"""
    ids_selected = get_selected_ids(domain, rank)
    output_filename = f"data/proteome-tree/{domain}-one_proteome_per_{rank}.trimmed.fa"
    with open(output_filename, 'w') as outfile:
        for acc in acc2seq_trimmed:
            if acc in ids_selected:
                outfile.write(f">{acc}\n{acc2seq_trimmed[acc]}\n")
        for acc in acc2seq_seeds_trimmed:
            outfile.write(f">{acc}\n{acc2seq_seeds_trimmed[acc]}\n")

# write_hmmalign_fasta('fungal', 'order')
# write_hmmalign_fasta('fungal', 'family')
write_hmmalign_fasta('all', 'class')
# write_hmmalign_fasta('all', 'order')

# write_trimmed_fasta('fungal', 'order')
# write_trimmed_fasta('fungal', 'family')
write_trimmed_fasta('all', 'class')
# write_trimmed_fasta('all', 'order')