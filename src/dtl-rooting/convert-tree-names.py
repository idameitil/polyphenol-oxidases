import json
import pandas as pd

def make_protein2proteome_dict(json):
    protein2proteome = {}
    for entry in json:
        protein_acc = entry['metadata']['accession']
        if len(entry['proteome_subset']) > 1:
            print(entry)
        proteome_acc = entry['proteome_subset'][0]['accession']
        protein2proteome[protein_acc] = proteome_acc.upper()
    return protein2proteome

export_json = json.load(open('data/proteome-tree/export.json'))
protein2proteome = make_protein2proteome_dict(export_json)

selected_genomes_file = 'data/proteome-tree/selected_genomes.xlsx'
df = pd.read_excel(selected_genomes_file)
proteome2species = {}
for index, row in df.iterrows():
    proteome2species[row.genome_id] = row.species

# Make replace_strings
replace_strings = {protein: proteome2species[protein2proteome[protein]].replace(' ', '_') for protein in protein2proteome if protein2proteome[protein] in proteome2species}

tree_filename = 'data/proteome-tree/raxml/T8.raxml.bestTree'
with open(tree_filename, "r+") as infile:
    content = infile.readlines()
    new_content = []

    for line in content:
        new_line = line

        for word in replace_strings.items():
            new_line = new_line.replace(str(word[0]), str(word[1]))
        new_content.append(new_line)

    with open('new.nwk', "w") as outfile:
        for line in new_content:
            outfile.write(line)