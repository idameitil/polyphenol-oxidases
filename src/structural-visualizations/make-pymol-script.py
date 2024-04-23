import json

def show_conserved_string(conserved_positions, object_name):
    string = ""
    for position in conserved_positions:
        string += f'label n. CA and resi {position} and {object_name}, "%s-%s" % (resn, resi)\n'
    temp_string = f"select cons_{object_name}, "
    for position in conserved_positions:
        temp_string += f"resi {position} and {object_name} or "
    string += temp_string[:-4] + '\n'
    string += f"show licorice, cons_{object_name}"
    return string

def load_model_string(object_name, pdb):
    return f"\nload {pdb}, {object_name}\n"

def pymol_object_name(acc):
    return f"{acc}"

# Load conserved residue
conserved_residues_filename = 'data/compare-architectures/conserved_0.95.json'
with open(conserved_residues_filename) as infile:
    conserved_residues = json.load(infile)

object_names = []

script = "@src/pymol-visualization/nicify.pml\n"

entries = [
    {'acc': '2p3x', 'family': 'a_plants', 'descriptive_name': 'PPOVv', 'domain_start_structure': 77, 'domain_end': 284, 'threshold':0.9}, 
    # {'acc': 'Q63JI5', 'family': 'b_cnidaria', 'descriptive_name': 'Q63JI5', 'domain_start_structure': 49, 'domain_end': 334, 'threshold':0.9}, 
    {'acc': '2y9x', 'family': 'c_long_fungal', 'descriptive_name': 'AbPPO3', 'domain_start_structure': 52, 'domain_end': 308, 'threshold':0.95}, 
    {'acc': '1wx2', 'family': 'd_bacteria', 'descriptive_name': 'TyrSc', 'domain_start_structure': 29, 'domain_end': 226, 'threshold':0.95}, 
    # {'acc': '1wx2', 'family': 'k_bacteria2', 'descriptive_name': 'TyrSc', 'domain_start_structure': 29, 'domain_end': 226, 'threshold':0.95}, 
    {'acc': '5m8l', 'family': 'e_chordata', 'descriptive_name': 'TyrHs', 'domain_start_structure': 184, 'domain_end': 416, 'threshold':0.95},
    {'acc': 'V3ZAB2', 'family': 'f_mollusc', 'descriptive_name': 'V3ZAB2', 'domain_start_structure': 130, 'domain_end': 300, 'threshold':0.93},
    {'acc': 'D0N318', 'family': 'h_oomycota', 'descriptive_name': 'D0N318', 'domain_start_structure': 77, 'domain_end': 274, 'threshold':0.99},
    {'acc': '4j3p', 'family': 'i_short_fungal', 'descriptive_name': 'AoCO4', 'domain_start_structure': 93, 'domain_end': 322, 'threshold':0.88},
    {'acc': 'F4PFF7', 'family': 'j_zoopagomycota', 'descriptive_name': 'F4PFF7', 'domain_start_structure': 68, 'domain_end': 247, 'threshold':0.92},
]

for entry in entries:
    model_path = f"data/compare-architectures/pdbs/{entry['acc']}.pdb"
    script += load_model_string(entry['acc'], model_path)
    script += show_conserved_string(conserved_residues[f"{entry['acc']} {entry['family']}"], entry['acc'])
    object_names.append(entry['acc'])

### Remove other chains ###
script += "\nremove chain B or chain C or chain D or chain E or chain F or chain G or chain H or chain I or chain J\n"

### ALIGN ###
script += f"""
cealign 2p3x, 2y9x
cealign 2p3x, 1wx2
cealign 2p3x, 5m8l
cealign 2p3x, V3ZAB2
cealign 2p3x, D0N318
cealign 2p3x, 4j3p
cealign 2p3x, D0N318
cealign 2p3x, F4PFF7
"""

# Nicify
script += "@src/structural-visualizations/nicify.pml\n"

output_filename = f"src/structural-visualizations/conserved-residues.pml"
with open (output_filename, 'w') as outfile:
    outfile.write(script)