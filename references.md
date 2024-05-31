# Removed genomes
Aureococcus anophagefferens and Ectocarpus siliculosus were removed. Almost all of the hits had low e-values and or short lengths, and they had no phylum.

Adiantum capillus-veneris was removed. Hits are very short.

Nitrososphaeraceae archaeon, de to hits har for lav e-værdi.

# Mount HPC
sshfs idamei@login1.gbar.dtu.dk:/work3/idamei hpc

sudo diskutil umount force hpc

# PDB sequences
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

# Tmux
List sessions: `tmux list-sessions`
Kill session: `tmux kill-session -t mysession`
Resume session: `tmux attach -d -t mrbayes`.

Check if the gpu is available: `nvtop -d 1`

#
Blåmusling mange tyrosinaser, adhesiveness

# Proteome analysis
Sequences from proteomes: 13381
Unique species with proteomes: 2472
Unique orders with proteomes: 335
Unique fungal orders with proteomes: 335
Fungal sequences from proteomes: 0
Unique fungal families with proteomes: 0
Sequences from selected fungal species: 654
Sequences from selected all species: 2254

# Albert
ssh dtu -tt ssh hal
Avogadro: adding hydrogens, molecular manipulations (also chemcraft)

grimme (xtb), optimization
also used for getting charges of Cu and Os
crest

orca, BP86

# Discussion with Caio 31 Jan.
Applications of oxydized phenols: Enhance the antioxidant value of lignin. Change the propoerties. Bacteria can only open rings of diphenols.

# hydroxybenzoic acid (H unit)

# Short fungal
Nikolaivits 2018 characterised XP_003666010.1
Gasparetti 2009 characterised 4J3P  (crystal is in another paper) (there are several later characterizations of this one)
Oates 2021 activity on lignin dimer CAI4219580.1
Frommhagen 2017 MtPPO7 XP_003662515.1

# CAG28310.1
Was deleted because it was super long. It seems to be several proteins.

# Installation of Interproscan on the HPC
See `Interproscan_installation.md`.

# Using efetch
`efetch -db protein -format fasta -input [accession_list] > [fasta_filename]`

# Phylogenetic trees
Try RAxML (as described in GH16 paper)

# Nitrosation
hydroxyannilinase = O-amino phenoloxidase

# Seeds
c2092_g1_i1 is from “A multi-omics approach to lignocellulolytic enzyme discovery reveals a new ligninase activity from Parascedosporium putredinis NO1”
c2092_g1_i1 has a blast hit with 100% coverage and identity except for one gap: CAI4219580.1

sp|C7FF04|PPO3_AGABI very well-studied, has a structure
sp|C7FF04|PPO3_AGABI has a blast hit with 100% identity and coverage: C7FF04.1

# Names
type-3 copper proteins and coupled binuclear cupper (CBC) proteins are not completely the same thing.

Laccases have a type-3 copper binding domain, but also a type-2 and type 1 domain

# The non-catalytic domains
- “The latent enzyme can be activated by limited proteolysis with proteinase K which cleaves the polypeptide chain after K382” (Pretzler 2017)
- “Latent and active abPPO4 mushroom tyrosinase cocrystallized with hexatungstotellurate(VI) in a single crystal” (Mauracher 2014)
- “The latent state of tyrosinase (pro-tyrosinase) is composed of the central domain and the C-terminal domain. The latter domain blocks the entrance into the active site by means of a “placeholder” residue which penetrates the active site similarly to the substrate or an inhibitor”, Kanteev, 2015
- “The N-terminal domain is a transit peptide which determines the final location of the enzyme and then undergoes proteolytic removal”