import os

dir = "data/mrbayes/0816/clades/members"
clades = os.listdir(dir)
outdir = "data/mrbayes/0816/clades/alignments-nofragments" 
for clade in clades:
    os.system(f'linsi --thread 7 data/mrbayes/0816/clades/fastas-nofragments/{clade}.fa > {outdir}/{clade}-linsi.fa')
