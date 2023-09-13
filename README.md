# Data collection

Caio sent a fasta file with the sequences for the proteins he has been working on: `data/caios-PPOs.fa`.

Caio sent a file with the PDB structures mentioned in Kanteev it al, 2015: `data/PPOs-with-a-3D-structure.xlsx`. The PDB IDs were extracted into a list: `data/pdb-ids`.

The sequences from the pdbs were downloaded at https://www.rcsb.org/downloads/fasta: `data/rcsb_pdb_20230911061612.fasta`.

The two fastas were combined into `data/pdbs-and-caios.fa`.

The seeds were cleaned in `data/pdbs-and-caios-cleaned.fa`. 
- 2Y9X_2 was removed, because it was called "LECTIN-LIKE FOLD PROTEIN". 2Y9W_1 and 2Y9W_2 were removed - they are the same sequence as 2Y9X. Only one is native state and the other is tropolone bound.
- 1WX2_2 was removed, because it was called "MelC".
- 3W6W was removed. The sequence is identical to 3W6Q. 6JU5, 6JU6 etc. were removed. They were almost identical to 3W6Q.
- 4J3P, 4J3Q, 5OR4, 5OR3, 4J3R and 6GSG were identical except for 4J3Q, which is missing the start. Only 4J3P was kept.

To be decided
- 4OUA and 5M6B: The two chains in 4OUA are identical at the N-terminus, but 4OUA_2 has an additional C-terminal domain. 5M6B is very similar to 4OUA_2.
- 3HHS: according to the paper, this PPO is a heterodimer of 2 homologous polypeptide chains PPO1 and PPO2 (Li, 2009). I guess we should keep both?

A seed table was created `data/seeds.tsv`.