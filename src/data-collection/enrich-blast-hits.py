import pandas as pd
from Bio import SeqIO

def read_interproscan(filename):
    df_interproscan = pd.read_csv(filename, sep='\t', names=['protein_accession', 'md5', 'sequence_length', \
                                                                'analysis', 'signature_accession', 'signature_discription',\
                                                                'start_location', 'stop_location', 'score', 'status',\
                                                                'date', 'interpro_annotations_accession', \
                                                                'interpro_annotations_description'])
    return df_interproscan

# Read interproscan
chunk2interpro_df = {}
unique_pfams_all = set()
N = 20
for i in range(0, N):
    df = read_interproscan(f"data/interproscan-blast-hits/chunk{'%02d' % i}.interproscan")
    chunk2interpro_df[i] = df
    unique_pfams = df[df.analysis == 'Pfam'].signature_accession.unique()
    unique_pfams_all.update(unique_pfams)

# Get pfam data
mydata = {}
count = 0
for i in range(0, N):
    df_interproscan = chunk2interpro_df[i]
    for acc in df_interproscan['protein_accession'].unique():
        count += 1
        print(count)
        mydata[acc] = []
        for pfam_name in unique_pfams_all:
            hits = df_interproscan[(df_interproscan.protein_accession == acc) & (df_interproscan.signature_accession == pfam_name)]
            if hits.empty:
                mydata[acc].append('')
            elif len(hits) == 1:
                mydata[acc].append(f"{hits.start_location.values[0]}-{hits.stop_location.values[0]}")
            else:
                string = ''
                for row, index in hits.iterrows():
                    string += f"{hits.start_location.values[0]}-{hits.stop_location.values[0]},"
                mydata[acc].append(string[:-1])

# Get sequences
count = 0
fasta_sequences = SeqIO.parse(open('data/blast/unique-hits-1e-15-length150-1000-cd-hit65.fasta'), 'fasta')
for fasta in fasta_sequences:
    if fasta.id in mydata:
        mydata[fasta.id].append(str(fasta.seq))
    else:
        count += 1
print(count)

column_names = list(unique_pfams_all)
column_names.append('sequence')
df_pfam_positions = pd.DataFrame.from_dict(mydata, orient='index', columns=column_names)
df_pfam_positions.index.name = 'id'
df_pfam_positions.to_csv("data/blast/unique-hits-1e-15-length150-1000-cd-hit65-enriched.tsv", sep='\t')