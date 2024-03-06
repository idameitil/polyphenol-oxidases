import json
import sys
sys.path.insert(0, 'src/data-collection')
from common import make_interproscan_df
import pandas as pd
import os

json_object = json.load(open('data/proteome-tree/export.json'))
acc2domain, unique_domains = {}, set()

for protein in json_object:
    protein_accession = protein['metadata']['accession']
    if os.path.isfile(f"data/pfam/api/{protein_accession}_data1.json"):
        with open(f"data/pfam/api/{protein_accession}_data1.json") as infile:
            data1 = json.load(infile)
        with open(f"data/pfam/api/{protein_accession}_data2.json") as infile:
            data2 = json.load(infile)
    for signature_description in data1:
        if data1[signature_description]['source_database'] not in ['phobius', 'signalp_g+', 'signalp_e']:
            continue
        signature_accession = data1[signature_description]['accession']
        domain_full_name = f"domain_{signature_accession}_{signature_description}"
        unique_domains.add(domain_full_name)
        start_location = data1[signature_description]['locations'][0]['fragments'][0]['start']
        stop_location = data1[signature_description]['locations'][0]['fragments'][0]['end']
        location_string = f"{start_location}-{stop_location}"
        if len(data1[signature_description]['locations']) > 1:
            for location in data1[signature_description]['locations'][1:]:
                next_start = location['fragments'][0]['start']
                next_stop = location['fragments'][0]['end']
                location_string += f",{next_start}-{next_stop}"
        if protein_accession not in acc2domain:
            acc2domain[protein_accession] = {domain_full_name: location_string}
        else:
            acc2domain[protein_accession][domain_full_name] = location_string
    # data2
    for result in data2['results']:
        if result['metadata']['source_database'] != 'pfam':
            continue
        signature_accession = result['metadata']['accession']
        signature_description = result['metadata']['name']
        if len(result['proteins']) > 1:
            print(result['proteins'])
        if result['proteins'][0]['entry_protein_locations'] is None:
            continue
        start_location = result['proteins'][0]['entry_protein_locations'][0]['fragments'][0]['start']
        stop_location = result['proteins'][0]['entry_protein_locations'][0]['fragments'][0]['end']
        domain_full_name = f"domain_{signature_accession}_{signature_description}"
        unique_domains.add(domain_full_name)
        location_string = f"{start_location}-{stop_location}"
        if len(result['proteins'][0]['entry_protein_locations']) > 1:
            for location in result['proteins'][0]['entry_protein_locations'][1:]:
                next_start = location['fragments'][0]['start']
                next_stop = location['fragments'][0]['end']
                location_string += f",{next_start}-{next_stop}"
        if protein_accession not in acc2domain:
            acc2domain[protein_accession] = {domain_full_name: location_string}
        else:
            acc2domain[protein_accession][domain_full_name] = location_string

df_interproscan = make_interproscan_df(acc2domain, unique_domains)

# Read enriched unique hits
df_unique_hits_enriched = pd.read_csv("data/pfam/protein-matching-PF00264.tsv", sep = '\t')

# Merge
merged_df = pd.merge(df_unique_hits_enriched, df_interproscan, how='left', on='protein_accession')
merged_df.to_csv('data/pfam/protein-matching-PF00264-interproscan2.tsv', sep='\t', index=False)