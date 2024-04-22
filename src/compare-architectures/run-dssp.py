from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
import os, sys

p = PDBParser()

def get_SS_string(dssp_object, discountinous):
    SS_string = 'type,end\n'
    previous = ''
    position = 0
    for residue in dssp_object:
        position += 1
        SS_translation = {'H': 'h', 'B': 'l', 'E': 's', 'G': 'h', 'I': 'h', 'T': 'l', 'S': 'l', '-': 'l'}
        simplified_secondary_structure = SS_translation[residue[2]]
        if position in discountinous:
            position += (discountinous[position] - position + 1)
        if simplified_secondary_structure == previous:
            continue
        if previous == 'h':
            SS_string += f"h,{position-1}\n"
        elif previous == 's':
            SS_string += f"s,{position-1}\n"
        elif previous == 'l':
            SS_string += f"l,{position-1}\n"
        previous = simplified_secondary_structure
    return SS_string

def write_ss_string(id, model_filename, outfilename):
    structure = p.get_structure(id, model_filename)
    model = structure[0]
    dssp = DSSP(model, model_filename)
    string = get_SS_string(dssp, wanted[id])
    with open(outfilename, 'w') as outfile:
        outfile.write(string)

# Values are discontinous intervals
wanted = {'4j3p': {1:4, 38:45}, 
          '5m8l': {1:25}, 
          '1wx2': {}, 
          '2y9x': {}, 
          '2p3x': {},
          'F4PFF7': {},
          'V3ZAB2': {},
          'D0N318': {},
          'Q63JI5': {},
          }
for id in wanted:
    print(id)
    write_ss_string(id, f"data/compare-architectures/pdbs/{id}.pdb", f"data/compare-architectures/architecture-tables-dssp/{id}.csv")
