from Bio import SeqIO
import json
import pandas as pd
import sys

def make_json_dict(json_list):
    json_dict = {}
    for entry in json_list:
        json_dict[entry['metadata']['accession']] = entry
    return json_dict

def make_protein2proteome_dict(json):
    protein2proteome = {}
    for entry in json:
        protein_acc = entry['metadata']['accession']
        if len(entry['proteome_subset']) > 1:
            print(entry)
        proteome_acc = entry['proteome_subset'][0]['accession']
        protein2proteome[protein_acc] = proteome_acc.upper()
    return protein2proteome

class proteomeData():

    df = pd.read_csv('data/proteome-tree/proteome-data.tsv', sep='\t')
    json_list = json.load(open('data/pfam/protein-matching-PF00264.json'))
    json_dict = make_json_dict(json_list)
    export_json = json.load(open('data/proteome-tree/export.json'))
    protein2proteome = make_protein2proteome_dict(export_json)

    def __init__(self, domain, rank):
        self.domain = domain
        self.rank = rank
        self.selected_proteomes2species = self.get_selected_proteomes()
        self.accs_from_selected_proteomes = self.get_accs_from_selected_proteomes()
        self.filtered_sequences = self.filter_sequences()

    def get_selected_proteomes(self): 
        if self.domain == 'all' and self.rank == 'class':
            df = pd.read_excel('data/proteome-tree/class-representatives.xlsx')
        elif self.domain == 'fungi' and self.rank == 'order':
            df = pd.read_excel('data/proteome-tree/fungal-order-representatives.xlsx')
        selected_proteomes2species = {} 
        for index, row in df.iterrows():
            selected_proteomes2species[row.genome_id] = row.species
        return selected_proteomes2species

    def get_accs_from_selected_proteomes(self):
        accs_from_selected_proteomes = list()
        for protein_acc in self.protein2proteome:
            if self.protein2proteome[protein_acc] in self.selected_proteomes2species:
                accs_from_selected_proteomes.append(protein_acc)
        return accs_from_selected_proteomes
    
    def filter_sequences(self):
        fasta_sequences = SeqIO.parse('data/pfam/protein-matching-PF00264-shortheaders.fasta', 'fasta')
        filtered_sequences = {}
        count = 0
        count_succeded = 0
        for fasta in fasta_sequences:
            if fasta.id not in self.accs_from_selected_proteomes:
                continue
            count += 1
            for entry_location in self.json_dict[fasta.id]['entries'][0]['entry_protein_locations']:
                score = entry_location['score']
                start = entry_location['fragments'][0]['start']
                end = entry_location['fragments'][0]['end']
                if score > 1e-20:
                    continue
                filtered_sequences[fasta.id] = f"{start}-{end}"
                count_succeded += 1
        print(f"Total: {count}")
        print(f"Succeeded: {count_succeded}")
        return filtered_sequences
    
    def write_selected_sequences(self):
        filename = f'data/proteome-tree/selected-sequences-{self.domain}-{self.rank}.txt'
        with open(filename, 'w') as outfile:
            for id in self.filtered_sequences:
                outfile.write(f"{id},{self.filtered_sequences[id]}\n")

all_class = proteomeData(domain='all', rank='class')
all_class.write_selected_sequences()