from Bio import SeqIO

fasta = "data/proteome-tree/sequences-all-class-filtered-andseeds.trimmed.fa"

def get_ids(fastas):
    ids = {}
    for fasta in fastas:
        fasta_sequences = SeqIO.parse(fasta, "fasta")
        for fasta_sequence in fasta_sequences:
            acc = fasta_sequence.id.split("/")[0]
            if acc in ids:
                if fasta_sequence.id.replace("/", ".") not in ids[acc]:
                    ids[acc].append(fasta_sequence.id.replace("/", "."))
            else:
                ids[acc] = [fasta_sequence.id.replace("/", ".")]
    return ids


ids = get_ids([fasta])


def translate(infile_name, outfile_name, sep):
    with open(infile_name) as infile, open(outfile_name, "w") as outfile:
        flag = False
        for line in infile:
            if not flag:
                outfile.write(line)
                if line.strip() == "DATA":
                    flag = True
            else:
                acc = line.split(sep)[0]
                if acc in ids:
                    for entry in ids[acc]:
                        outfile.write(f"{entry}{sep}{sep.join(line.split(sep)[1:])}")


itol_files = [
    {"name": "domain-combined", "sep": ","},
    {"name": "domain-combined-1500", "sep": ","},
    {"name": "uniprot-kingdom-strip", "sep": "\t"},
    {"name": "uniprot-phylum-strip", "sep": "\t"},
    {"name": "uniprot-class-text", "sep": ","},
    {"name": "uniprot-species-text", "sep": ","},
    {"name": "uniprot-adapted-strip", "sep": "\t"},
]

for file in itol_files:
    name = file["name"]
    infile = f"data/itol-label-files/{name}.txt"
    outfile = f"data/itol-label-files/{name}-location.txt"
    translate(infile, outfile, sep=file["sep"])
