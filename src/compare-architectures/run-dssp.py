from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
import os, sys

p = PDBParser()

def get_SS_string(dssp_object, discountinous, id):
    SS_string = 'type,end\n'
    previous = ''
    position = 0
    count = 0
    for residue in dssp_object:
        if id == '1lnl':
            count += 1
            if count < 4:
                continue
        print(residue)
        previous_position = position
        position += 1
        SS_translation = {'H': 'h', 'B': 'l', 'E': 's', 'G': 'h', 'I': 'h', 'T': 'l', 'S': 'l', '-': 'l'}
        if position in discountinous:
            simplified_secondary_structure = 'u'
            position = discountinous[position] + 1
        else:
            simplified_secondary_structure = SS_translation[residue[2]]
        if simplified_secondary_structure == previous:
            continue
        if previous == 'h':
            SS_string += f"h,{previous_position}\n"
        elif previous == 's':
            SS_string += f"s,{previous_position}\n"
        elif previous == 'l':
            SS_string += f"l,{previous_position}\n"
        elif previous == 'u':
            SS_string += f"u,{previous_position}\n"
        previous = simplified_secondary_structure
    return SS_string

def write_ss_string(id, model_filename, outfilename):
    structure = p.get_structure(id, model_filename)
    model = structure[0]
    dssp = DSSP(model, model_filename)
    string = get_SS_string(dssp, wanted[id], id)
    with open(outfilename, 'w') as outfile:
        outfile.write(string)

# Values are discontinous intervals
wanted = {'4j3p': {1:4, 38:45}, 
          '5m8l': {1:24}, 
          '1wx2': {1:1}, 
          '2y9x': {1:1}, 
          '2p3x': {},
          'F4PFF7': {},
          'V3ZAB2': {},
          'D0N318': {},
          'Q63JI5': {},
          '1js8': {1:2502, 2637:2644},
          'A0A7M6DMJ9': {},
          '1lnl': {},
          '6z1s': {1:14, 49:52, 192:195, 338:340, 406:425},
          'A0A137PAT2': {},
          '5zrd': {1:3},
          'O76708': {}
          }
for id in wanted:
    print(id)
    write_ss_string(id, f"data/compare-architectures/pdbs/{id}.pdb", f"data/compare-architectures/architecture-tables-dssp/{id}.csv")
