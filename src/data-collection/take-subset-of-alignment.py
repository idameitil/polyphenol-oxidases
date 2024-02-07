from Bio import SeqIO

# Read filtered
filtered_fasta_filename = 'data/pfam/protein-matching-PF00264-shortheaders-filtered.fasta'
fasta_sequences = SeqIO.parse(filtered_fasta_filename, 'fasta')
accessions_include = []
for fasta in fasta_sequences:
    accessions_include.append(fasta.id)

# Read fungi
fungi_fasta_filename = 'data/pfam/protein-matching-PF00264-fungi-shortheaders.fasta'
fasta_sequences = SeqIO.parse(fungi_fasta_filename, 'fasta')
accessions_fungi = []
for fasta in fasta_sequences:
    accessions_fungi.append(fasta.id)

# Read alignment
alignment_filename = 'data/pfam/PF00264.alignment.uniprot-cleaned.fa'
fasta_sequences = SeqIO.parse(alignment_filename, 'fasta')

# Write output
output_filename_all = 'data/pfam/PF00264.alignment.uniprot-cleaned-filtered.fa'
output_filename_all_withoutgaps = 'data/pfam/PF00264.alignment.uniprot-cleaned-filtered-withoutgaps.fa'
output_filename_fungi = 'data/pfam/PF00264.alignment.uniprot-cleaned-filtered-fungi.fa'
output_filename_fungi_withoutgaps = 'data/pfam/PF00264.alignment.uniprot-cleaned-filtered-fungi-withoutgaps.fa'

accessions_done = []
with open(output_filename_all, 'w') as outfile_all, \
open(output_filename_all_withoutgaps, 'w') as outfile_all_withoutgaps, \
open(output_filename_fungi, 'w') as outfile_fungi, \
open(output_filename_fungi_withoutgaps, 'w') as outfile_fungi_withoutgaps:
    for fasta in fasta_sequences:
        if fasta.id.split('.')[0] in accessions_include:
            outfile_all.write(f">{fasta.id}\n{fasta.seq}\n")
            outfile_all_withoutgaps.write(f">{fasta.id}\n{fasta.seq.replace('-', '')}\n")
            # Write fungal fasta
            if fasta.id.split('.')[0] in accessions_fungi:
                outfile_fungi.write(f">{fasta.id}\n{fasta.seq}\n")
                outfile_fungi_withoutgaps.write(f">{fasta.id}\n{fasta.seq.replace('-', '')}\n")