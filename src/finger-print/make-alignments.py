import os

dir = "data/mrbayes/all/clades/members"
clades = os.listdir(dir)
outdir = "data/mrbayes/all/clades/alignments" 
for clade in clades:
    os.system(f'linsi --thread 7 data/mrbayes/all/clades/fastas/{clade}.fa > {outdir}/{clade}-linsi.fa')