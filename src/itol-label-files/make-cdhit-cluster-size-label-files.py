from cdhit_reader import read_cdhit
from colour import Color

def make_itol_file(infile, outfile):
    representative2size = {}
    sizes = []
    for cluster in read_cdhit(infile):
        representative_name = cluster.refname.split('.')[0]
        size = len(cluster.sequences)
        sizes.append(size)
        representative2size[representative_name] = size
    min_size = min(sizes)
    max_size = max(sizes)
    n = max_size - min_size
    colors = list(Color("white").range_to(Color("black"),n+1))
    with open(outfile, 'w') as outfile:
        header = header = f"DATASET_COLORSTRIP\nSEPARATOR TAB\nDATASET_LABEL\tcluster_size\nCOLOR\t#ff0000\nDATA\n"
        outfile.write(header)
        for representative_name in representative2size:
            size = representative2size[representative_name]
            color = colors[size-min_size].hex
            outfile.write(f"{representative_name}\t{color}\t{size}\n")
    
# Fungal
print('*Fungal 0.5*\n')
fungi_5_infile = "data/pfam/PF00264.alignment.uniprot-cleaned-filtered25-fungi-withoutgaps-cdhit0.5.fasta.clstr"
fungi_5_outfile = "data/pfam/PF00264.alignment.uniprot-cleaned-filtered25-fungi-withoutgaps-cdhit0.5.fasta-size-itol.txt"
make_itol_file(fungi_5_infile, fungi_5_outfile)

print('*Fungal 0.6*\n')
fungi_6_infile = "data/pfam/PF00264.alignment.uniprot-cleaned-filtered25-fungi-withoutgaps-cdhit0.6.fasta.clstr"
fungi_6_outfile = "data/pfam/PF00264.alignment.uniprot-cleaned-filtered25-fungi-withoutgaps-cdhit0.6.fasta-size-itol.txt"
make_itol_file(fungi_6_infile, fungi_6_outfile)

# All kingdoms
print('*All kingdoms 0.5*\n')
all_5_infile = "data/pfam/PF00264.alignment.uniprot-cleaned-filtered25-withoutgaps-cdhit0.5.fasta.clstr"
all_5_outfile = "data/pfam/PF00264.alignment.uniprot-cleaned-filtered25-withoutgaps-cdhit0.5.fasta-size-itol.txt"
make_itol_file(all_5_infile, all_5_outfile)

print('*All kingdoms 0.4*\n')
all_4_infile = "data/pfam/PF00264.alignment.uniprot-cleaned-filtered25-withoutgaps-cdhit0.4.fasta.clstr"
all_4_outfile = "data/pfam/PF00264.alignment.uniprot-cleaned-filtered25-withoutgaps-cdhit0.4.fasta-size-itol.txt"
make_itol_file(all_4_infile, all_4_outfile)