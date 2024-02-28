from Bio import SeqIO

fasta_selected = SeqIO.parse('data/proteome-tree/fungal-one_proteome_per_order.fa', 'fasta')
ids_selected = [entry.id for entry in fasta_selected]

#########################################################
# For HMMalign tree
#########################################################

fasta_alignment = SeqIO.parse('data/pfam/PF00264.alignment.uniprot-cleaned.fa', 'fasta')

output_filename = "data/proteome-tree/fungal-one_proteome_per_order.hmmalign.fa"
with open(output_filename, 'w') as outfile:
    for fasta in fasta_alignment:
        if fasta.id.split('.')[0] in ids_selected:
            outfile.write(f">{fasta.id}\n{fasta.seq}\n")
        
#########################################################
# For mafft tree (make trimmed fasta)
#########################################################
            
fasta_alignment = SeqIO.parse('data/pfam/PF00264.alignment.uniprot-nogaps.fa', 'fasta')

output_filename = "data/proteome-tree/fungal-one_proteome_per_order.trimmed.fa"
with open(output_filename, 'w') as outfile:
    for fasta in fasta_alignment:
        if fasta.id.split('.')[0] in ids_selected:
            outfile.write(f">{fasta.id}\n{fasta.seq}\n")