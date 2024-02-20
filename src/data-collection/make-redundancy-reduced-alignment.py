from Bio import SeqIO

fasta_alignment = SeqIO.parse('data/pfam/PF00264.alignment.uniprot-cleaned.fa', 'fasta')

# fastas_all_40 = SeqIO.parse('data/pfam/PF00264.alignment.uniprot-cleaned-filtered-withoutgaps-cdhit0.4.fasta', 'fasta')
# fastas_all_50 = SeqIO.parse('data/pfam/PF00264.alignment.uniprot-cleaned-filtered-withoutgaps-cdhit0.5.fasta', 'fasta')
# fastas_fungi_50 = SeqIO.parse('data/pfam/PF00264.alignment.uniprot-cleaned-filtered-fungi-withoutgaps-cdhit0.5.fasta', 'fasta')
# fastas_fungi_60 = SeqIO.parse('data/pfam/PF00264.alignment.uniprot-cleaned-filtered-fungi-withoutgaps-cdhit0.6.fasta', 'fasta')
fastas_all_40 = SeqIO.parse('data/pfam/PF00264.alignment.uniprot-cleaned-filtered25-withoutgaps-cdhit0.4.fasta', 'fasta')
fastas_all_50 = SeqIO.parse('data/pfam/PF00264.alignment.uniprot-cleaned-filtered25-withoutgaps-cdhit0.5.fasta', 'fasta')
fastas_fungi_50 = SeqIO.parse('data/pfam/PF00264.alignment.uniprot-cleaned-filtered25-fungi-withoutgaps-cdhit0.5.fasta', 'fasta')
fastas_fungi_60 = SeqIO.parse('data/pfam/PF00264.alignment.uniprot-cleaned-filtered25-fungi-withoutgaps-cdhit0.6.fasta', 'fasta')

ids_all_40 = [entry.id for entry in fastas_all_40]
ids_all_50 = [entry.id for entry in fastas_all_50]
ids_fungi_50 = [entry.id for entry in fastas_fungi_50]
ids_fungi_60 = [entry.id for entry in fastas_fungi_60]

# output_filename_all_40 = "data/pfam/PF00264.alignment.uniprot-cleaned-filtered-cdhit0.4.fasta"
# output_filename_all_50 = "data/pfam/PF00264.alignment.uniprot-cleaned-filtered-cdhit0.5.fasta"
# output_filename_fungi_50 = "data/pfam/PF00264.alignment.uniprot-cleaned-filtered-fungi-cdhit0.5.fasta"
# output_filename_fungi_60 = "data/pfam/PF00264.alignment.uniprot-cleaned-filtered-fungi-cdhit0.6.fasta"
output_filename_all_40 = "data/pfam/PF00264.alignment.uniprot-cleaned-filtered25-cdhit0.4.fasta"
output_filename_all_50 = "data/pfam/PF00264.alignment.uniprot-cleaned-filtered25-cdhit0.5.fasta"
output_filename_fungi_50 = "data/pfam/PF00264.alignment.uniprot-cleaned-filtered25-fungi-cdhit0.5.fasta"
output_filename_fungi_60 = "data/pfam/PF00264.alignment.uniprot-cleaned-filtered25-fungi-cdhit0.6.fasta"
with open(output_filename_all_40, 'w') as outfile_all_40,\
open(output_filename_all_50, 'w') as outfile_all_50,\
open(output_filename_fungi_50, 'w') as outfile_fungi_50,\
open(output_filename_fungi_60, 'w') as outfile_fungi_60:
    for fasta in fasta_alignment:
        if fasta.id in ids_all_40:
            outfile_all_40.write(f">{fasta.id}\n{fasta.seq}\n")
        if fasta.id in ids_all_50:
            outfile_all_50.write(f">{fasta.id}\n{fasta.seq}\n")
        if fasta.id in ids_fungi_50:
            outfile_fungi_50.write(f">{fasta.id}\n{fasta.seq}\n")
        if fasta.id in ids_fungi_60:
            outfile_fungi_60.write(f">{fasta.id}\n{fasta.seq}\n")