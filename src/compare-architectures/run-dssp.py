from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
import os, sys

p = PDBParser()

def get_SS_string(dssp_object):
    SS_string = 'type,end\n'
    previous = ''
    for residue in dssp_object:
        SS_translation = {'H': 'h', 'B': 'l', 'E': 's', 'G': 'h', 'I': 'h', 'T': 'l', 'S': 'l', '-': 'l'}
        simplified_secondary_structure = SS_translation[residue[2]]
        position = residue[0]
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
    string = get_SS_string(dssp)
    with open(outfilename, 'w') as outfile:
        outfile.write(string)

wanted = ['4j3p', '5m8l', '1wx2', '2y9x', '2p3x']
for id in wanted:
    write_ss_string(id, f"data/compare-architectures/pdbs/{id}.pdb", f"data/compare-architectures/architecture-tables-dssp/{id}.csv")
