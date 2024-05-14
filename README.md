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
`/work3/idamei/bin/my_interproscan/interproscan-5.64-96.0/interproscan.sh -appl Pfam,SignalP_EUK,SignalP_GRAM_NEGATIVE,SignalP_GRAM_POSITIVE,Phobius -i /work3/idamei/polyphenol-oxidases/seeds.fa -f tsv -o /work3/idamei/polyphenol-oxidases/seeds.interproscan`

Download the output file:
`scp idamei@transfer.gbar.dtu.dk:/work3/idamei/polyphenol-oxidases/seeds.interproscan data/`

## Make fasta with only pfam domain
A fasta file with only the PF00264 and PF00372 domains of the seeds was made by running: `python src/interproscan/make-pfam-domain-fasta.py`. This creates the file `data/seeds-pfam-domains.fa`.

## Enriching seeds
The seed table is enriched with taxonomy by running `python src/data-collection/enrich-seeds.py`. This creates the file `data/seeds-enriched.tsv`.

# Phylogenetic tree
A phylogenetic tree is made on the mafft webserver and saved in `data/seeds.ph`.

# iTOL label files
To make iTOL label files, run `python src/itol-label-files/make-itol-label-files.py`

# Structural visualizations
To show the structures run `pymol src/structural-visualizations/show-models.pml`.

# Downloading Pfam family
Fasta of the Pfam entries was downloaded from `https://www.ebi.ac.uk/interpro/entry/pfam/PF00264/protein/UniProt/#table` and saved in `data/pfam/protein-matching-PF00264.fasta`.

A json file with metadata was also downloaded and saved in `data/pfam/protein-matching-PF00264.json`.

(Not used) Aligned sequences were downloaded `data/pfam/PF00264.alignment.uniprot`.

The pfam HMM was downloaded at `https://www.ebi.ac.uk/interpro/entry/pfam/PF00264/curation/` and saved in `data/pfam/PF00264.hmm`.

A json file was downloaded with all proteome data (e-mail correspondence uniprot): `data/proteome-tree/export.json`. To download the domain data through the interproscan api, run `python src/proteome-tree/download-domain-data-api.py`. This downloads the files to `data/pfam/api`.

# Processing of uniprot fasta
A fasta file with only the uniprot ID in the header was generated using find and replace regex (find ">([A-Z0-9]+).+", replace with ">$1") on `data/pfam/protein-matching-PF00264.fasta`: `data/pfam/protein-matching-PF00264-shortheaders.fasta`




# Retrieve taxonomy
To make a tsv file with taxonomy for each uniprot hit, run: `python src/data-collection/make-taxonomy-file.py`. This creates the file `data/pfam/protein-matching-PF00264-fungi.tsv`.

### Enrich fungal uniprot hits with domain architecture
To enrich the uniprot hits with pfam data, run `python src/data-collection/enrich-uniprot-hits-interproscan.py`. This creates the file `data/pfam/protein-matching-PF00264-interproscan.tsv`. (note that only the reduced hits are enriched, so there are many lines without pfam info)

New: `python src/proteome-tree/enrich-proteome-hits-interproscan.py`. This creates the file `data/pfam/protein-matching-PF00264-interproscan2.tsv`.

# Eggnog
Eggnog was run at `http://eggnog-mapper.embl.de/` with this input file: `data/pfam/PF00264.alignment.uniprot-cleaned-filtered-withoutgaps.fa`.

The output files were saved in: `data/eggnog`.

A table with only the OGs is made by running: `python src/eggnog/parse-eggnog.py`.

# Proteome tree

## Get proteome data
A json file with the proteomes matching PF00264 was downloaded from pfam: `data/pfam/proteome-matching-PF00264.json`.

A json file with proteome metadata was downloaded from uniprot `https://www.uniprot.org/proteomes?query=*`: `data/proteome-tree/proteomes_AND_proteome_type_1_2024_02_28.json`.

The above to files are combined to create the table by running: `python src/proteome-tree/make-proteome-table.py`. This creates the file `data/proteome-tree/proteome-data.tsv`.

## Manually selected genomes
(Old)
A file with manually selected genomes was saved in `data/proteome-tree/selected_genomes.xlsx`

A file with excluded species was saved in `data/proteome-tree/exclude`.
(New)
Files with manually selected genomes were saved in `data/proteome-tree/class-representatives.xlsx` and `data/proteome-tree/fungal-order-representatives.xlsx`

## Select proteomes and get sequences
In order to select the proteomes to include and get the PPO sequences from those proteomes, run `python src/proteome-tree/get-proteome-sequences.py`. This creates the set of files `data/proteome-tree/selected-proteomes-ids-fungi-order.txt` (containing the taxids for the included proteomes) and `data/proteome-tree/fungal-one_proteome_per_order.fa` (containing the PPO sequences for the selected proteomes) for each set of parameters.

## Get aligned selected sequences
In order to make fasta files with the aligned selected sequences, run `python get-aligned-selected-sequences.py`. This creates the files `data/proteome-tree/all-one_proteome_per_class.hmmalign.fa` and data/proteome-tree/all-one_proteome_per_class.trimmed.fa etc.

## Make trees
To make hmmalign tree, run `raxml-ng --msa data/proteome-tree/fungal-one_proteome_per_order.hmmalign.fa --model JTT+G4 --prefix data/proteome-tree/raxml/T1 --threads 7 --seed 2 --blopt nr_safe`

To make hmmalign tree, run `raxml-ng --msa data/proteome-tree/all-one_proteome_per_class.hmmalign.fa --model JTT+G4 --prefix data/proteome-tree/raxml/T2 --threads 7 --seed 2 --blopt nr_safe`

To make hmmalign tree, run `raxml-ng --msa data/proteome-tree/fungal-one_proteome_per_family.hmmalign.fa --model JTT+G4 --prefix data/proteome-tree/raxml/T3 --threads 7 --seed 2 --blopt nr_safe`

To make hmmalign tree, run `raxml-ng --msa data/proteome-tree/all-one_proteome_per_class.hmmalign.fa --model JTT+G4 --prefix data/proteome-tree/raxml/T4 --threads 7 --seed 2 --blopt nr_safe`

With even shorter: `raxml-ng --msa data/proteome-tree/fungal-one_proteome_per_order.hmmalign-evenshorter.fa --model JTT+G4 --prefix data/proteome-tree/raxml/T5 --threads 7 --seed 2`

With even shorter: `raxml-ng --msa data/proteome-tree/fungal-one_proteome_per_family.hmmalign-evenshorter.fa --model JTT+G4 --prefix data/proteome-tree/raxml/T6 --threads 7 --seed 2`

### Bootstrap
`raxml-ng --msa data/proteome-tree/all-one_proteome_per_class.hmmalign.fa --model JTT+G4 --prefix data/proteome-tree/raxml/T7 --threads 7 --seed 2 --all --bs-trees 100`

## With mafft
### All, one per class
Run mafft: `linsi --thread 7 data/proteome-tree/all-one_proteome_per_class.trimmed.fa > data/proteome-tree/all-one_proteome_per_class.trimmed.linsi.fa`.

Remove columns with more than 90% gaps: `seqconverter -I fasta -O fasta --remfracgapcols 0.9 data/proteome-tree/all-one_proteome_per_class.trimmed.linsi.fa > data/proteome-tree/all-one_proteome_per_class.trimmed.linsi-0.9.fa`

Convert to nexus: `seqconverter -I fasta -O nexus data/proteome-tree/all-one_proteome_per_class.trimmed.linsi-0.9.fa > data/proteome-tree/all-one_proteome_per_class.trimmed.linsi-0.9.nexus`

raxml-ng --msa data/proteome-tree/all-one_proteome_per_class.trimmed.linsi-0.9.fa --model JTT+G4 --prefix data/proteome-tree/raxml/T13 --threads 7 --seed 2 --all --bs-trees 200

### Fungal, one per order
Run mafft: `linsi --thread 7 data/proteome-tree/fungal-one_proteome_per_order.trimmed.fa > data/proteome-tree/fungal-one_proteome_per_order.trimmed.linsi.fa`

Remove columns with more than 90% gaps: `seqconverter -I fasta -O fasta --remfracgapcols 0.9 data/proteome-tree/fungal-one_proteome_per_order.trimmed.linsi.fa > data/proteome-tree/fungal-one_proteome_per_order.trimmed.linsi-0.9.fa`

Convert to nexus: `seqconverter -I fasta -O nexus data/proteome-tree/fungal-one_proteome_per_order.trimmed.linsi-0.9.fa > data/proteome-tree/fungal-one_proteome_per_order.trimmed.linsi-0.9.nexus`

## DTL rooting
An input file was made with the species tree and gene tree (in the gene tree, the genes are named by species): `data/dtl/input-dtl` by running `src/dtl-rooting/convert-tree-names.py`.

ranger-dtl was run with varying parameters, eg: `../ranger-dtl/CorePrograms/OptRoot.mac -i data/dtl/input-dtl -o data/dtl/dtl128.out -D 2 -L 1 -T 8`

# Lignin degraders
The list of lignin degrading species was saved in `data/lignin-degraders`.

# Clades
Files with proteins in each clade (determined from iTOL) were saved in `data/mrbayes/all/clades`.

A csv file is made by running `python src/proteome-tree/clade-annotations.py`.

iTOL heatmap file is made by running `python src/itol-label-files/make-itol-label-files.py` and is shown on the species tree.

# MrBayes
On hal, go to the folder and start tmux: `tmux new -s mrbayes`

Start the job: `mpirun -np 6 mb gpu_run.nexus > log.txt`

Detach: control+b followed by d

(it seems like you don't have to do this step) qTo continue a mrbayes run, copy everything from the .ckp~ file into the bottom of the aignment nexus file. In the run nexus file, change the number of generations to the total number of generations you want and add append=yes, fx: `mcmc append=yes ngen=60000000 samplefreq=1000 nchains=3 file=out.nex`.

The run with all kingdoms, one per order was saved in `data/mrbayes/all`.

# Tmux
List sessions: `tmux list-sessions`
Kill session: `tmux kill-session -t mysession`
Resume session: `tmux attach -d -t mrbayes`.

Check if the gpu is available: `nvtop -d 1`

# Fingerprint
The ids in each clade were saved in `data/mrbayes/all/clades`.

Fasta files are generated by running `python src/finger-print/make-clade-fastas.py`.

# Pymol
To make pymol script with aligned structures and conserved residues, run `python src/structural-visualizations/make-pymol-script.py`.

To run the pymol script, run `pymol src/structural-visualizations/conserved-residues.pml`.