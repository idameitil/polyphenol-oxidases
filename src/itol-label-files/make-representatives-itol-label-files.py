from cdhit_reader import read_cdhit
import sys

def make_itol_file(infile, outfile, seed_list):
    representative2members = {}
    for cluster in read_cdhit(infile):
        for member in cluster.sequences:
            if member.name in seed_list:
                if cluster.refname == member.name and len(cluster) == 1:
                    print(f'sequence {member.name} not in tree\n')
                elif cluster.refname == member.name and len(cluster) > 1:
                    print(f'sequence {member.name} is the new representative. Check cluster')
                    print(f"{cluster.name} refSequence={cluster.refname} size={len(cluster)}")
                    print(f" {member.name} ({member.length}) identity={member.identity}% {'(Reference sequence)' if member.is_ref else ''}\n")
                elif cluster.refname in seed_list:
                    print(f"sequence {member.name} is represented by another seed {cluster.refname}\n")
                else:
                    representative_name = cluster.refname.split('.')[0]
                    if representative_name in representative2members:
                        representative2members[representative_name].append(member.name)
                    else:
                        representative2members[representative_name] = [member.name]
    with open(outfile, 'w') as outfile:
        header = "DATASET_TEXT\nSEPARATOR COMMA\nDATASET_LABEL,representative-sequences\nCOLOR,#000000\nDATA\n"
        outfile.write(header)
        for representative in representative2members:
            members_string = ';'.join(representative2members[representative])
            outfile.write(f"{representative},{members_string},-1,#000000,bold,10,0\n")

# Fungal
fungal_seed_names = ['abPPO4', 'AoMelB', 'AoCO4','CgPPO1473','CgPPO266','MtPPO809','TtPPO010','PpPPOc2029','MtPPO7','AbPPO3']

print('*Fungal 0.5*\n')
fungi_5_infile = "data/pfam/PF00264.alignment.uniprot-cleaned-filtered-fungi-withoutgaps-withseeds-cdhit0.5.fasta.clstr"
fungi_5_outfile = "data/pfam/PF00264.alignment.uniprot-cleaned-filtered-fungi-withoutgaps-withseeds-cdhit0.5.fasta-itol.txt"
make_itol_file(fungi_5_infile, fungi_5_outfile, fungal_seed_names)

print('*Fungal 0.6*\n')
fungi_6_infile = "data/pfam/PF00264.alignment.uniprot-cleaned-filtered-fungi-withoutgaps-withseeds-cdhit0.6.fasta.clstr"
fungi_6_outfile = "data/pfam/PF00264.alignment.uniprot-cleaned-filtered-fungi-withoutgaps-withseeds-cdhit0.6.fasta-itol.txt"
make_itol_file(fungi_6_infile, fungi_6_outfile, fungal_seed_names)

# All kingdoms
all_seed_names = ['ScTyr','BmTyr','VvPPO','IbCO','abPPO4','AoMelB','proPOMc_1','proPOMc_2','ProPOc','Lobster_Hc','Crab_Hc',
              'EdHc','RvHc','AoCO4','CgPPO1473','CgPPO266','MtPPO809','TtPPO010','PpPPOc2029','MtPPO7',
              'AbPPO3','NpsF','SlPPO1','SlPPO2','MdPPO1','MdPPO2','MdPPO3','CgAUS1','TyrHs']

print('*All kingdoms 0.5*\n')
all_5_infile = "data/pfam/PF00264.alignment.uniprot-cleaned-filtered-withoutgaps-withseeds-cdhit0.5.fasta.clstr"
all_5_outfile = "data/pfam/PF00264.alignment.uniprot-cleaned-filtered-withoutgaps-withseeds-cdhit0.5.fasta-itol.txt"
make_itol_file(all_5_infile, all_5_outfile, all_seed_names)

print('*All kingdoms 0.4*\n')
all_4_infile = "data/pfam/PF00264.alignment.uniprot-cleaned-filtered-withoutgaps-withseeds-cdhit0.4.fasta.clstr"
all_4_outfile = "data/pfam/PF00264.alignment.uniprot-cleaned-filtered-withoutgaps-withseeds-cdhit0.4.fasta-itol.txt"
make_itol_file(all_4_infile, all_4_outfile, all_seed_names)

# Filtering e-25
print('*Fungal e-25 0.5*\n')
fungi_25_5_infile = "data/pfam/PF00264.alignment.uniprot-cleaned-filtered25-fungi-withoutgaps-withseeds-cdhit0.5.fasta.clstr"
fungi_25_5_outfile = "data/pfam/PF00264.alignment.uniprot-cleaned-filtered25-fungi-withoutgaps-withseeds-cdhit0.5.fasta-itol.txt"
make_itol_file(fungi_25_5_infile, fungi_25_5_outfile, fungal_seed_names)

print('*Fungal e-25 0.6*\n')
fungi_25_6_infile = "data/pfam/PF00264.alignment.uniprot-cleaned-filtered25-fungi-withoutgaps-withseeds-cdhit0.6.fasta.clstr"
fungi_25_6_outfile = "data/pfam/PF00264.alignment.uniprot-cleaned-filtered25-fungi-withoutgaps-withseeds-cdhit0.6.fasta-itol.txt"
make_itol_file(fungi_25_6_infile, fungi_25_6_outfile, fungal_seed_names)

print('*All kingdoms e-25 0.5*\n')
all_25_5_infile = "data/pfam/PF00264.alignment.uniprot-cleaned-filtered25-withoutgaps-withseeds-cdhit0.5.fasta.clstr"
all_25_5_outfile = "data/pfam/PF00264.alignment.uniprot-cleaned-filtered25-withoutgaps-withseeds-cdhit0.5.fasta-itol.txt"
make_itol_file(all_25_5_infile, all_25_5_outfile, all_seed_names)

print('*All kingdoms e-25 0.4*\n')
all_25_4_infile = "data/pfam/PF00264.alignment.uniprot-cleaned-filtered25-withoutgaps-withseeds-cdhit0.4.fasta.clstr"
all_25_4_outfile = "data/pfam/PF00264.alignment.uniprot-cleaned-filtered25-withoutgaps-withseeds-cdhit0.4.fasta-itol.txt"
make_itol_file(all_25_4_infile, all_25_4_outfile, all_seed_names)