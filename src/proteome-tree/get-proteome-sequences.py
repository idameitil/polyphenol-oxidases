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
    df_manually_selected = pd.read_excel('data/proteome-tree/selected_genomes.xlsx')
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
        for fasta in fasta_sequences:
            if fasta.id not in self.accs_from_selected_proteomes:
                continue
            if len(self.json_dict[fasta.id]['entries']) > 1:
                # print(f"{fasta.id} has several entries")
                continue
            if len(self.json_dict[fasta.id]['entries'][0]['entry_protein_locations']) > 1:
                # print(f"{fasta.id} has several entry_protein_locations")
                continue
            if len(self.json_dict[fasta.id]['entries'][0]['entry_protein_locations'][0]['fragments']) > 1:
                # print(f"{fasta.id} has several fragments")
                continue
            hit_length = self.json_dict[fasta.id]['entries'][0]['entry_protein_locations'][0]['fragments'][0]['end'] \
                - self.json_dict[fasta.id]['entries'][0]['entry_protein_locations'][0]['fragments'][0]['start']
            score = self.json_dict[fasta.id]['entries'][0]['entry_protein_locations'][0]['score']
            if len(fasta.seq) > 100 and len(fasta.seq) < 1000 and score < 1e-20 and hit_length > 100:
                filtered_sequences[fasta.id] = fasta.seq
        return filtered_sequences
    
    def fasta_filename(self):
        if self.domain == 'fungi':
            return f'data/proteome-tree/fungal-one_proteome_per_{self.rank}.fa'
        elif self.domain == 'all':
            return f'data/proteome-tree/all-one_proteome_per_{self.rank}.fa'

    def write_fasta(self):
        with open(self.fasta_filename(), 'w') as outfile:
            for acc in self.filtered_sequences:
                proteome_id = self.protein2proteome[acc]
                species = self.selected_proteomes2species[proteome_id]
                # outfile.write(f">{species}_{acc}\n{self.filtered_sequences[acc]}\n")
                outfile.write(f">{acc}\n{self.filtered_sequences[acc]}\n")

    def write_proteomes_txt(self):
        filename = f'data/proteome-tree/selected-proteomes-ids-{self.domain}-{self.rank}.txt'
        with open(filename, 'w') as outfile:
            for id in self.selected_proteomes2species:
                outfile.write(f"{id}\n")

fungi_order = proteomeData(domain='fungi', rank='order')
fungi_order.write_fasta()
fungi_order.write_proteomes_txt()

all_class = proteomeData(domain='all', rank='class')
all_class.write_fasta()
all_class.write_proteomes_txt()
