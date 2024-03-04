from Bio import SeqIO
import json
import pandas as pd
import sys

def make_json_dict(json_list):
    json_dict = {}
    for entry in json_list:
        json_dict[entry['metadata']['accession']] = entry
    return json_dict

class proteomeData():

    df = pd.read_csv('data/proteome-tree/proteome-data.tsv', sep='\t')
    json_list = json.load(open('data/pfam/protein-matching-PF00264.json'))
    json_dict = make_json_dict(json_list)

    def __init__(self, domain, rank):
        self.domain = domain
        self.rank = rank
        self.selected_proteome_taxids = self.get_selected_proteome_taxids()
        self.accs_from_selected_proteomes = self.get_accs_from_selected_proteomes()
        self.filtered_sequences = self.filter_sequences()

    def get_selected_proteome_taxids(self):
        if self.domain == 'fungi':
            unique_taxons = self.df[self.df['kingdom'] == 'Fungi'][self.rank].unique()
        elif self.domain == 'all':
            unique_taxons = self.df[self.rank].unique()
        selected_proteome_taxids = []
        for taxon in unique_taxons:
            if pd.isna(taxon):
                continue
            best_proteome = self.df[self.df[self.rank] == taxon]['score'].idxmax()
            best_proteome_taxid = self.df.iloc[[best_proteome]].taxid.item()
            selected_proteome_taxids.append(best_proteome_taxid)
        return selected_proteome_taxids

    def get_accs_from_selected_proteomes(self):
        accs_from_selected_proteomes = list()
        for entry in self.json_list:
            if int(entry['metadata']['source_organism']['taxId']) in self.selected_proteome_taxids:
                accs_from_selected_proteomes.append(entry['metadata']['accession'])
        return accs_from_selected_proteomes

    def filter_sequences(self):
        fasta_sequences = SeqIO.parse('data/pfam/protein-matching-PF00264-shortheaders.fasta', 'fasta')
        filtered_sequences = {}
        for fasta in fasta_sequences:
            if fasta.id not in self.accs_from_selected_proteomes:
                continue
            if len(self.json_dict[fasta.id]['entries']) > 1:
                print(f"{fasta.id} has several entries")
                continue
            if len(self.json_dict[fasta.id]['entries'][0]['entry_protein_locations']) > 1:
                print(f"{fasta.id} has several entry_protein_locations")
                continue
            if len(self.json_dict[fasta.id]['entries'][0]['entry_protein_locations'][0]['fragments']) > 1:
                print(f"{fasta.id} has several fragments")
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
                outfile.write(f">{acc}\n{self.filtered_sequences[acc]}\n")

    def write_proteomes_txt(self):
        filename = f'data/proteome-tree/selected-proteomes-taxids-{self.domain}-{self.rank}.txt'
        with open(filename, 'w') as outfile:
            for taxid in self.selected_proteome_taxids:
                outfile.write(f"{taxid}\n")

fungi_order = proteomeData(domain='fungi', rank='order')
fungi_order.write_fasta()
fungi_order.write_proteomes_txt()

fungi_family = proteomeData(domain='fungi', rank='family')
fungi_family.write_fasta()
fungi_family.write_proteomes_txt()

all_order = proteomeData(domain='all', rank='order')
all_order.write_fasta()
all_order.write_proteomes_txt()

all_class = proteomeData(domain='all', rank='class')
all_class.write_fasta()
all_class.write_proteomes_txt()
