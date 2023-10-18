import pandas as pd

# Read enriched unique hits
df_unique_hits_enriched = pd.read_csv("data/blast/unique-hits-enriched.tsv", sep = '\t')

# Read interproscan output
N = 20
unique_pfams = set()
acc2pfam = {}
for i in range(0, N):
    filename = f"data/interproscan-blast-hits/chunk{'%02d' % i}.interproscan"
    with open(filename, 'r') as infile:
        for line in infile:
            splitted_line = line.split('\t')
            protein_accession, md5, sequence_length, analysis, signature_accession, \
            signature_description, start_location, stop_location, score, status,\
            date, interpro_annotations_accession, interpro_annotations_description \
            = splitted_line
            if analysis == 'Pfam':
                pfam_full_name = f"{signature_accession}_{signature_description}"
                unique_pfams.add(pfam_full_name)
                if protein_accession not in acc2pfam:
                    acc2pfam[protein_accession] = {pfam_full_name: f"{start_location}-{stop_location}"}
                else:
                    if pfam_full_name not in acc2pfam[protein_accession]:
                        acc2pfam[protein_accession][pfam_full_name] = f"{start_location}-{stop_location}"
                    else:
                        acc2pfam[protein_accession][pfam_full_name] += f",{start_location}-{stop_location}"
                        print(acc2pfam[protein_accession][pfam_full_name])

# Make interproscan dataframe
mydata = {}
count = 0
for acc in acc2pfam:
    count += 1
    print(count)
    mydata[acc] = []
    for pfam_full_name in unique_pfams:
        if pfam_full_name not in acc2pfam[acc]:
            mydata[acc].append('')
        else:
            mydata[acc].append(acc2pfam[acc][pfam_full_name])
column_names = list(unique_pfams)
df_positions = pd.DataFrame.from_dict(mydata, orient='index', columns=column_names)
df_positions.index.name = 'protein_accession'

# Merge
merged_df = pd.merge(df_unique_hits_enriched, df_positions, how='right', on='protein_accession')
merged_df.to_csv('data/blast/unique-hits-1e-60-length150-1000-cd-hit65-enriched.tsv', sep='\t')