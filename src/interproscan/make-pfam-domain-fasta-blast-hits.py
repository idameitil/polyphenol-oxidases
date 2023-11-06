import pandas as pd

# Blast hits
enriched_hits_df = pd.read_csv("data/blast/unique-hits-1e-60-length150-1000-cd-hit65-enriched.tsv", sep='\t')

output_filename = 'data/blast/unique-hits-1e-60-length150-1000-cd-hit65-pfam-domains.fa'
with open(output_filename, 'w') as outfile:
    for index, row in enriched_hits_df.iterrows():
        if not pd.isnull(row["domain_PF00264_Common central domain of tyrosinase"]):
            if ',' in row["domain_PF00264_Common central domain of tyrosinase"]:
                continue
            else:
                start, stop = row["domain_PF00264_Common central domain of tyrosinase"].split('-')
                pfam_domain = row.seq[int(start)-1:int(stop)-1]
                outfile.write(f">{row.protein_accession}\n{pfam_domain}\n")