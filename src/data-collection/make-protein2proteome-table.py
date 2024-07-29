import json
import pandas as pd
from common import get_taxon

def make_protein_dict(json):
    protein2proteome = {}
    for entry in json:
        protein_acc = entry["metadata"]["accession"]
        if len(entry["proteome_subset"]) > 1:
            print(entry)
        proteome_acc = entry["proteome_subset"][0]["accession"]
        taxid = entry["metadata"]['source_organism']["taxId"]
        species = get_taxon(taxid=taxid, rank='species')
        protein2proteome[protein_acc] = {'proteome_acc': proteome_acc.upper(), 'species': species}
    return protein2proteome

export_json = json.load(open("data/proteome-tree/export.json"))
protein2proteome = make_protein_dict(export_json)

df = pd.DataFrame.from_dict(protein2proteome).T
df.index.name = 'protein_acc'

df.to_csv('data/proteome-tree/protein2proteome.tsv', sep='\t')