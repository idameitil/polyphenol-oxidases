import pandas as pd
import random


def make_pfam_label_file(df):
    outfilename = f"data/itol-label-files/blast-hits-pfam.txt"

    with open(outfilename, "w") as file:
        header = f"DATASET_DOMAINS\nSEPARATOR COMMA\nDATASET_LABEL,Domains\nCOLOR,#ff0000\nDATA\n"
        file.write(header)
        pfam2color = dict()
        count = 0
        for index, row in df.iterrows():
            try:
                file.write(f"{row.id},{len(row.sequence)}")
            except:
                count += 1
                continue
            for pfam_name in df.columns[1:-1]:
                if pfam_name not in pfam2color:
                    color = '#' + "%06x" % random.randint(0, 0xFFFFFF)
                    pfam2color[pfam_name] = color
                if not pd.isnull(row[pfam_name]):
                    # Several occurences of this pfam in the sequence
                    if ',' in row[pfam_name]:
                        occurences = row[pfam_name].split(',')
                        for occurence in occurences:
                            start, stop = occurence.split('-')
                            file.write(f",OC|{start}|{stop}|{pfam2color[pfam_name]}|{pfam_name}")
                    # One occurence of the pfam in the sequence
                    else:
                        start, stop = row[pfam_name].split('-')
                        file.write(f",OC|{start}|{stop}|{pfam2color[pfam_name]}|{pfam_name}")
            file.write('\n')
        print(count)

df = pd.read_csv('data/blast/unique-hits-1e-15-length150-1000-cd-hit65-enriched.tsv', sep='\t')
make_pfam_label_file(df)