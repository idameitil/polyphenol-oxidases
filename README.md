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

## Fasta file
A fasta file is created by running `python src/data-collection/make-fasta.py`. This creates the file `data/seeds.fa`.

## Enriching seeds
The seed table is enriched with taxonomy by running `python src/data-collection/enrich-seeds.py`. This creates the file `data/seeds-enriched.tsv`.

# MSA
An MSA is made by running `mafft data/seeds.fa > data/seeds-mafft.fa`

# Phylogenetic tree
A phylogenetic tree is made on the mafft webserver and saved in `data/seeds.ph`.

# iTOL label files
To make the iTOL label files, run `python src/make-itol-label-files.py`

# Structural visualizations
To show the structures run `pymol src/structural-visualizations/show-models.pml`.