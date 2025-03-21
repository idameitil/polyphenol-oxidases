from Bio import SeqIO
import json
import pandas as pd
import sys


def make_json_dict(json_list):
    json_dict = {}
    for entry in json_list:
        json_dict[entry["metadata"]["accession"]] = entry
    return json_dict


def make_protein2proteome_dict(json):
    protein2proteome = {}
    for entry in json:
        protein_acc = entry["metadata"]["accession"]
        if len(entry["proteome_subset"]) > 1:
            print(entry)
        proteome_acc = entry["proteome_subset"][0]["accession"]
        protein2proteome[protein_acc] = proteome_acc.upper()
    return protein2proteome


class proteomeData:

    df = pd.read_csv("data/proteome-tree/proteome-data.tsv", sep="\t")
    json_list = json.load(open("data/pfam/protein-matching-PF00264.json"))
    json_dict = make_json_dict(json_list)
    export_json = json.load(open("data/proteome-tree/export.json"))
    protein2proteome = make_protein2proteome_dict(export_json)

    def __init__(self, domain, rank, score_threshold):
        self.domain = domain
        self.rank = rank
        self.score_threshold = score_threshold
        self.selected_proteomes2species = self.get_selected_proteomes()
        self.accs_from_selected_proteomes = self.get_accs_from_selected_proteomes()
        self.selected_sequences, self.filtered_out_sequences = self.filter_sequences(self.score_threshold)

    def get_selected_proteomes(self):
        if self.domain == "all" and self.rank == "class":
            df = pd.read_excel("data/proteome-tree/class-representatives.xlsx")
        elif self.domain == "fungi" and self.rank == "order":
            df = pd.read_excel("data/proteome-tree/fungal-order-representatives.xlsx")
        elif self.domain == "fungi-with-lignin-degraders" and self.rank == "order":
            df = pd.read_excel("data/proteome-tree/fungal-order-representatives-andlignindegraders.xlsx")
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

    def filter_sequences(self, score_threshold):
        fasta_sequences = SeqIO.parse(
            "data/pfam/protein-matching-PF00264-shortheaders.fasta", "fasta"
        )
        selected_sequences = []
        filtered_out_sequences = []
        count = 0
        for fasta in fasta_sequences:
            if fasta.id not in self.accs_from_selected_proteomes:
                continue
            count += 1
            for entry_location in self.json_dict[fasta.id]["entries"][0]["entry_protein_locations"]:
                score = entry_location["score"]
                start = entry_location["fragments"][0]["start"]
                end = entry_location["fragments"][0]["end"]
                if score > score_threshold:
                    filtered_out_sequences.append(
                    {"acc": fasta.id, "position": f"{start}-{end}"}
                )
                else:
                    selected_sequences.append(
                        {"acc": fasta.id, "position": f"{start}-{end}"}
                    )
        print(f"Total: {count}")
        print(f"Succeeded: {len(selected_sequences)}")
        return selected_sequences, filtered_out_sequences

    def write_SequenceIds(self, outfilename, included, excluded=[]):
        with open(outfilename, 'w') as outfile:
            for entry in included:
                if entry not in excluded:
                    outfile.write(f"{entry['acc']},{entry['position']}\n")

    def write_all_sequences(self):
        filename = f"data/proteome-tree/sequenceIds-{self.domain}-{self.rank}-notFiltered.txt"
        included = self.selected_sequences + self.filtered_out_sequences 
        self.write_SequenceIds(filename, included)

    def write_selected_sequences(self):
        filename = f"data/proteome-tree/sequenceIds-{self.domain}-{self.rank}-filtered.txt"
        self.write_SequenceIds(filename, self.selected_sequences)

    def write_filtered_out_sequences(self):
        filename = f"data/proteome-tree/sequenceIds-{self.domain}-{self.rank}-removed.txt"
        self.write_SequenceIds(filename, self.filtered_out_sequences)

    def write_selected_sequences_nooverlap(self, excluded):
        filename = f"data/proteome-tree/sequenceIds-{self.domain}-{self.rank}-filtered-noOverlap.txt"
        self.write_SequenceIds(filename, self.selected_sequences, excluded)

    def write_filtered_out_sequences_nooverlap(self, excluded):
        filename = f"data/proteome-tree/sequenceIds-{self.domain}-{self.rank}-removed-noOverlap.txt"
        self.write_SequenceIds(filename, self.filtered_out_sequences, excluded)

all_class = proteomeData(domain="all", rank="class", score_threshold=1e-20)
all_class.write_selected_sequences()
all_class.write_filtered_out_sequences()

all_fungal = proteomeData(domain="fungi", rank="order", score_threshold=1e-20)
all_fungal.write_selected_sequences()
all_fungal.write_filtered_out_sequences()

all_fungal.write_filtered_out_sequences_nooverlap(excluded=all_class.filtered_out_sequences)
all_fungal.write_selected_sequences_nooverlap(excluded=all_class.selected_sequences)
all_fungal.write_all_sequences()

fungal_with_lignin = proteomeData(domain="fungi-with-lignin-degraders", rank="order", score_threshold=1e-20)
fungal_with_lignin.write_filtered_out_sequences_nooverlap(excluded=all_class.filtered_out_sequences)
fungal_with_lignin.write_selected_sequences_nooverlap(excluded=all_class.selected_sequences)