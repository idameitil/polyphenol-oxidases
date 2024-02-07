# Login directly to hal
ssh dtu -tt ssh hal

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