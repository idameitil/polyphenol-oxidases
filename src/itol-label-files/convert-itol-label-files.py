from Bio import SeqIO

fasta = 'data/proteome-tree/all-one_proteome_per_class.trimmed.fa'

def get_ids(fasta):
    fasta_sequences = SeqIO.parse(fasta, 'fasta')
    ids = {}
    for fasta_sequence in fasta_sequences:
        acc = fasta_sequence.id.split('.')[0]
        if acc in ids:
            ids[acc].append(fasta_sequence.id.replace('/', '0'))
        else:
            ids[acc] = [fasta_sequence.id.replace('/', '0')]
    return ids

ids = get_ids(fasta)

def translate(infile_name, outfile_name, sep):
    with open(infile_name) as infile, open(outfile_name, 'w') as outfile:
        flag = False
        for line in infile:
            if not flag:
                outfile.write(line)
                if line.strip() == 'DATA':
                    flag = True
            else:
                acc = line.split(sep)[0]
                if acc in ids:
                    for entry in ids[acc]:
                        outfile.write(f"{entry}{sep}{sep.join(line.split(sep)[1:])}")

itol_files = [{'name': 'domain-combined', 'sep': ','}, 
              {'name': 'uniprot-kingdom-strip', 'sep': '\t'},
              {'name': 'uniprot-phylum-strip', 'sep': '\t'},
              {'name': 'uniprot-class-text', 'sep': ','},
              {'name': 'uniprot-species-text', 'sep': ','},
              {'name': 'uniprot-adapted-strip', 'sep': '\t'}
              ]

for file in itol_files:
    name = file['name']
    infile = f'data/itol-label-files/{name}.txt' 
    outfile = f'data/itol-label-files/{name}-location.txt' 
    translate(infile, outfile, sep = file['sep'])
