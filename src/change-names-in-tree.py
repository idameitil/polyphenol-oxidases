import re

with open('data/pdbs-and-caios-cleaned.ph', 'r') as file:
    tree = file.read()

tree_cleaned = re.sub(r'_\d\|Chain', '', tree)
with open('data/pdbs-and-caios-cleaned-correct-headers.ph', 'w') as outfile:
    outfile.write(tree_cleaned)