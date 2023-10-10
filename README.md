# Data collection

Caio sent a fasta file with the sequences for the proteins he has been working on.

Caio sent a file with the PDB structures mentioned in Kanteev it al, 2015: `data/PPOs-with-a-3D-structure.xlsx`.

The sequences from the pdbs were downloaded at https://www.rcsb.org/downloads/fasta.

The pdb sequences were cleaned:
- 2Y9X_2 was removed, because it was called "LECTIN-LIKE FOLD PROTEIN". 2Y9W_1 and 2Y9W_2 were removed - they are the same sequence as 2Y9X. Only one is native state and the other is tropolone bound. UPDATED: Both 2Y9X were removed because they are identical to sp|C7FF04|PPO3_AGABI.
- 1WX2_2 was removed, because it was called "MelC".
- 3W6W was removed. The sequence is identical to 3W6Q. 6JU5, 6JU6 etc. were removed. They were almost identical to 3W6Q.
- 4J3P, 4J3Q, 5OR4, 5OR3, 4J3R and 6GSG were identical except for 4J3Q, which is missing the start. Only 4J3P was kept.
- 4OUA and 5M6B: The two chains in 4OUA are identical at the N-terminus, but 4OUA_2 has an additional C-terminal domain. 5M6B is very similar to 4OUA_2. 4OUA_1 was kept.
- 6Z1S is identical to XP_003666010.1, so it was deleted

To be decided:
- 3HHS: according to the paper, this PPO is a heterodimer of 2 homologous polypeptide chains PPO1 and PPO2 (Li, 2009). I guess we should keep both?

## Seed table
A seed table was created manually: `data/seeds.tsv`.

## Aguilera data
The supplementary table with the dataset from the Aguilera data was downloaded and saved in `data/Aguilera-data/12862_2013_2332_MOESM1_ESM.xlsx`.

The sequences for the ones with accession numbers were retrieved with efetch and put in the table `data/Aguilera-data/aguilera_with_seq.xlsx`

The data was include in a table with the same setup as `data/seeds.tsv` - `data/Aguilera-data/aguilera-dataset.xlsx`.

## Fasta file
A fasta file is created by running `python src/data-collection/make-fasta.py`. This creates the file `data/seeds.fa`.

## Running Interproscan on seeds
Transfer the seeds file to the HPC: `scp data/seeds.fa idamei@transfer.gbar.dtu.dk:/work3/idamei/polyphenol-oxidases/`

On the HPC, run Interproscan:
`/work3/idamei/bin/my_interproscan/interproscan-5.64-96.0/interproscan.sh -i /work3/idamei/polyphenol-oxidases/seeds.fa -f tsv -o /work3/idamei/polyphenol-oxidases/seeds.interproscan`

Download the output file:
`scp idamei@transfer.gbar.dtu.dk:/work3/idamei/polyphenol-oxidases/seeds.interproscan data/`

## Make fasta with only pfam domain
A fasta file with only the PF00264 and PF00372 domains of the seeds was made by running: `python src/interproscan/make-pfam-domain-fasta.py`. This creates the file `data/seeds-pfam-domains.fa`.

## Enriching seeds
The seed table is enriched with taxonomy by running `python src/data-collection/enrich-seeds.py`. This creates the file `data/seeds-enriched.tsv`.

## Expanding with Blast
Copy fasta file to HPC: `scp data/seeds-pfam-domains.fa idamei@transfer.gbar.dtu.dk:/work3/idamei/polyphenol-oxidases/blast/polyphenol-oxidases.fasta`

On the HPC, run: `python3 /work3/idamei/src/blast.py polyphenol-oxidases`

And then: `sh /work3/idamei/polyphenol-oxidases/blast/submit.sh`

(OBS: maybe only download the new runs) When all jobs have finished, download the data: `scp -r idamei@transfer.gbar.dtu.dk:/work3/idamei/polyphenol-oxidases/blast data/`

### Parsing Blast results
To parse the blast expansion output files into json format, run `python src/data-collection/run-blastfilter.py`. This will create the files `data/blast/run/*/blast.js`.

### Make unique-hits file
To create `unique-hits.tsv`, run `python src/data-collection/make-unique-hits-file.py`. This creates the file `data/blast/unique-hits.tsv`.

# MSA
An MSA is made by running `mafft data/seeds.fa > data/seeds-mafft.fa`

# Phylogenetic tree
A phylogenetic tree is made on the mafft webserver and saved in `data/seeds.ph`.

# iTOL label files
To make the iTOL label files, run `python src/itol-label-files/make-itol-label-files.py`

# Structural visualizations
To show the structures run `pymol src/structural-visualizations/show-models.pml`.


