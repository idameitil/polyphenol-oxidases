import sys, os, json
import pandas as pd

df = pd.read_csv(f"data/seeds.tsv", sep='\t')

blast_path = f"data/blast/run"

accessions = [folder for folder in os.listdir(blast_path) if not folder.startswith('.')]

hits_best_pct = {}
for accession in accessions:
    if not accession in df.protein_accession.to_list():
        continue
    json_filename = f"{blast_path}/{accession}/blast.js"
    with open(json_filename, 'r') as infile:
        blast_hits = []
        for line in infile:
            if line.startswith('#'):
                continue
            blast_hit = json.loads(line.strip())
            acc = blast_hit['name2']
            pct = blast_hit['pct']
            evalue = blast_hit['evalue']
            if acc not in hits_best_pct:
                hits_best_pct[acc] = {'pct': pct, 'evalue': evalue}
            else:
                if pct > hits_best_pct[acc]['pct']:
                    hits_best_pct[acc] = {'pct': pct, 'evalue': evalue}

# Write file with unique hits and their best e-value
output_filename = f"data/blast/unique-hits.tsv"
pd.DataFrame.from_dict(hits_best_pct).T.to_csv(output_filename, sep = '\t', header=False)