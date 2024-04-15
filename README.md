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

## Expanding with Blast
Copy fasta file to HPC: `scp data/seeds-pfam-domains.fa idamei@transfer.gbar.dtu.dk:/work3/idamei/polyphenol-oxidases/blast/polyphenol-oxidases.fasta`

On the HPC, run: `python3 /work3/idamei/src/blast.py polyphenol-oxidases`

And then: `sh /work3/idamei/polyphenol-oxidases/blast/submit.sh`

(OBS: maybe only download the new runs) When all jobs have finished, download the data: `scp -r idamei@transfer.gbar.dtu.dk:/work3/idamei/polyphenol-oxidases/blast data/`

### Parsing Blast results
To parse the blast expansion output files into json format, run `python src/data-collection/run-blastfilter.py`. This will create the files `data/blast/run/*/blast.js`.

### Make unique-hits file
To create `unique-hits.tsv`, run `python src/data-collection/make-unique-hits-file.py`. This creates the file `data/blast/unique-hits.tsv`.

### Retrieving sequences, taxid for Blast hits
To retrieve the sequence and taxid for the blast hits, run:
`scp data/blast/unique-hits.tsv idamei@transfer.gbar.dtu.dk:/work3/idamei/polyphenol-oxidases/blast/`

Then, on the HPC run:
`qrsh`

`blastdbcmd -db nr -entry_batch /work3/idamei/polyphenol-oxidases/blast/unique-hits.tsv > /work3/idamei/polyphenol-oxidases/blast/unique-hits.fasta`

`blastdbcmd -db nr -entry_batch /work3/idamei/polyphenol-oxidases/blast/unique-hits.tsv -outfmt "%a ,%L ,%T ,%t ,%s" > /work3/idamei/polyphenol-oxidases/blast/unique-hits.csv` 

Then locally run:
`scp idamei@transfer.gbar.dtu.dk:/work3/idamei/polyphenol-oxidases/blast/unique-hits.fasta data/blast/`

`scp idamei@transfer.gbar.dtu.dk:/work3/idamei/polyphenol-oxidases/blast/unique-hits.csv data/blast/`

`data/blast/unique-hits.fasta` is a fasta file with exactly the entries in `data/blast/unique-hits.tsv`. `data/blast/unique-hits.csv` is a csv file containing taxids along with other info. It also includes identical sequences and therefore has more lines than `data/blast/unique-hits.tsv`.

### Enrich unique hits with taxonomy
To enrich all unique blast hits with taxonomy, run `python src/data-collection/enrich-blast-hits.py`. This creates the file `data/blast/unique-hits-enriched.tsv`.

### Filter blast hits
To filter the blast hits on e-value and length, run `python src/data-collection/filter-blast-hits/filter-blast-hits.py 1e-60`. This creates the files `data/blast/unique-hits-1e-60.fasta` and `data/blast/unique-hits-1e-60-length150-1000.fasta`.

To make the 1e-30 fasta, run `python src/data-collection/filter-blast-hits/filter-blast-hits.py 1e-30`.

`cd-hit -i data/blast/unique-hits-1e-60-length150-1000.fasta -c 0.65 -o data/blast/unique-hits-1e-60-length150-1000-cd-hit65.fasta`.

### Run Interproscan on blast hits
Copy the fasta to the HPC: `scp data/blast/unique-hits-1e-60-length150-1000-cd-hit65.fasta idamei@transfer.gbar.dtu.dk:/work3/idamei/polyphenol-oxidases/interproscan-blast/`

On the HPC:
Chunk fasta: `chunkfasta -c 20 -d polyphenol-oxidases/interproscan-blast polyphenol-oxidases/interproscan-blast/unique-hits-1e-60-length150-1000-cd-hit65.fasta`

Run this:
`qrsh`
`/work3/idamei/bin/my_interproscan/interproscan-5.64-96.0/interproscan.sh -appl Pfam,SignalP_EUK,SignalP_GRAM_NEGATIVE,SignalP_GRAM_POSITIVE,Phobius -i /work3/idamei/polyphenol-oxidases/interproscan-blast/chunk00.fa -f tsv -o /work3/idamei/polyphenol-oxidases/interproscan-blast/chunk00.interproscan`

`/work3/idamei/bin/my_interproscan/interproscan-5.64-96.0/interproscan.sh -appl Pfam,SignalP_EUK,SignalP_GRAM_NEGATIVE,SignalP_GRAM_POSITIVE,Phobius -i /work3/idamei/polyphenol-oxidases/interproscan-blast/chunk01.fa -f tsv -o /work3/idamei/polyphenol-oxidases/interproscan-blast/chunk01.interproscan`
etc.

Download the results: `scp -r idamei@transfer.gbar.dtu.dk:/work3/idamei/polyphenol-oxidases/interproscan-blast data`

(not working yet, because it runs with the wrong version of Java, when submitted to the queue)
Using jobscripts:
On the HPC:
`python3 polyphenol-oxidases/interproscan-blast/runinterproscan.py`

### Enrich selected hits with Pfam
To enrich the filtered blast hits with pfam data, run `python src/data-collection/enrich-blast-hits-interproscan.py`. This creates the file `data/blast/unique-hits-enriched-interproscan.tsv`.

### Extract pfam domains from blast hits
To extract the pfam domains from the blast hits, run `python src/interproscan/make-pfam-domain-fasta-blast-hits.py`. This creates the file `data/blast/unique-hits-1e-60-length150-1000-cd-hit65-pfam-domains.fa`.

# MSA
An MSA is made by running `mafft data/seeds.fa > data/seeds-mafft.fa`

# Phylogenetic tree
A phylogenetic tree is made on the mafft webserver and saved in `data/seeds.ph`.

# iTOL label files
To make iTOL label files, run `python src/itol-label-files/make-itol-label-files.py`

# Structural visualizations
To show the structures run `pymol src/structural-visualizations/show-models.pml`.

# Downloading Pfam family (uniprot members etc.)
Fasta of the Pfam entries was downloaded from `https://www.ebi.ac.uk/interpro/entry/pfam/PF00264/protein/UniProt/#table` and saved in `data/pfam/protein-matching-PF00264.fasta`.
`cd-hit -i data/pfam/protein-matching-PF00264.fasta -c 0.4 -n 2 -o data/pfam/protein-matching-PF00264-cdhit0.4.fasta`

A fasta file with only the uniprot ID in the header was generated using find and replace regex (find ">([A-Z0-9]+).+", replace with ">$1"): `data/pfam/protein-matching-PF00264-shortheaders.fasta`

A json file with metadata was also downloaded and saved in `data/pfam/protein-matching-PF00264.json`.

Aligned sequences were downloaded `data/pfam/PF00264.alignment.uniprot`.

The alignment is cleaned from columns with insertions by running python `src/data-collection/clean-hmm-alignment.py`. This created the fasta file `data/pfam/PF00264.alignment.uniprot-cleaned.fa`.

A trimmed fasta was made: `data/pfam/PF00264.alignment.uniprot-nogaps.fa` by removing all gaps from the alignment.

The pfam HMM was downloaded at `https://www.ebi.ac.uk/interpro/entry/pfam/PF00264/curation/` and saved in `data/pfam/PF00264.hmm`.

# Make tree of whole family
CD-HIT: `cd-hit -i data/pfam/protein-matching-PF00264-shortheaders.fasta -c 0.4 -n 2 -o data/pfam/protein-matching-PF00264-shortheaders-cdhit0.4.fasta`

HMMalign: `hmmalign --trim data/pfam/PF00264.hmm data/pfam/protein-matching-PF00264-shortheaders-cdhit0.4.fasta > data/pfam/protein-matching-PF00264-shortheaders-cdhit0.4.fa.hmmaligntrim`.

Convert output file to fasta online and make upper case (http://sequenceconversion.bugaco.com/converter/biology/sequences/stockholm_to_fasta.php): `data/pfam/protein-matching-PF00264-shortheaders-cdhit0.4.fa.hmmaligntrim.fasta`

Make tree: `raxml-ng --search1 --msa data/pfam/protein-matching-PF00264-shortheaders-cdhit0.4.fa.hmmaligntrim.fasta --model JTT+G4 --prefix data/pfam/raxml/T2 --threads 8 --seed 2`

# Make tree of whole family with length filter
Filter fasta by length: `python src/data-collection/length-filter.py`. This generates the file `data/pfam/protein-matching-PF00264-shortheaders-70-1000.fasta`.

CD-HIT: `cd-hit -i data/pfam/protein-matching-PF00264-shortheaders-70-1000.fasta -c 0.4 -n 2 -o data/pfam/protein-matching-PF00264-shortheaders-70-1000-cdhit0.4.fasta`

HMMalign: `hmmalign --trim data/pfam/PF00264.hmm data/pfam/protein-matching-PF00264-shortheaders-70-1000-cdhit0.4.fasta > data/pfam/protein-matching-PF00264-shortheaders-70-1000-cdhit0.4.fa.hmmaligntrim`.

Convert output file to fasta online and make upper case (http://sequenceconversion.bugaco.com/converter/biology/sequences/stockholm_to_fasta.php): `data/pfam/protein-matching-PF00264-shortheaders-20-1000-cdhit0.4.fa.hmmaligntrim.fasta`

Make tree: `raxml-ng --search1 --msa data/pfam/protein-matching-PF00264-shortheaders-70-1000-cdhit0.4.fa.hmmaligntrim.fasta --model JTT+G4 --prefix data/pfam/raxml/T3 --threads 8 --seed 2`

# Make tree of whole family with length, score and coverage filter
Filter fasta `python src/data-collection/filter.py`. This generates the file `data/pfam/protein-matching-PF00264-shortheaders-filtered.fasta`.

CD-HIT: `cd-hit -i data/pfam/protein-matching-PF00264-shortheaders-filtered.fasta -c 0.4 -n 2 -o data/pfam/protein-matching-PF00264-shortheaders-filtered-cdhit0.4.fasta`

HMMalign: `hmmalign --trim data/pfam/PF00264.hmm data/pfam/protein-matching-PF00264-shortheaders-filtered-cdhit0.4.fasta > data/pfam/protein-matching-PF00264-shortheaders-filtered-cdhit0.4.fasta.hmmaligntrim`.

Convert output file to fasta online and make upper case (http://sequenceconversion.bugaco.com/converter/biology/sequences/stockholm_to_fasta.php): `data/pfam/protein-matching-PF00264-shortheaders-filtered-cdhit0.4.fasta.hmmaligntrim.fasta`

Make tree: `raxml-ng --search1 --msa data/pfam/protein-matching-PF00264-shortheaders-filtered-cdhit0.4.fasta.hmmaligntrim.fasta --model JTT+G4 --prefix data/pfam/raxml/T5 --threads 8 --seed 2`

Aclust tree with full length was made on the HPC and saved in `data/trees/aclust-uniprot-al-kingdoms-filtered.nwk`. This file was used as input: `data/pfam/protein-matching-PF00264-shortheaders-filtered-cdhit0.4.fasta`.

An aclust tree with trimmed sequences was made on the HPC and saved in `data/trees/aclust-uniprot-al-kingdoms-filtered.nwk`. This file was used as input: `data/pfam/protein-matching-PF00264-shortheaders-filtered-cdhit0.4.fasta.hmmaligntrim.fasta.aclustinput`.

## Make tree based on pfam alignment with filtered hits all kingdoms
To make an alignment file with only the filtered redundancy reduced hits, run `src/data-collection/take-subset-of-alignment.py`. This creates the file `data/pfam/PF00264.alignment.uniprot-cleaned-filtered.fa`.

To make the tree, run: `raxml-ng --msa data/pfam/PF00264.alignment.uniprot-cleaned-filtered.fa --model JTT+G4 --prefix data/pfam/raxml/T7 --threads 7 --seed 2 --blopt nr_safe`.

### With cd-hit on the alignment subset
To take the subset of the alignment that is in the fasta, run: `python src/data-collection/take-subset-of-alignment.py`. This creates the files `data/pfam/PF00264.alignment.uniprot-cleaned-filtered.fa`, `data/pfam/PF00264.alignment.uniprot-cleaned-filtered-withoutgaps.fa`, `data/pfam/PF00264.alignment.uniprot-cleaned-filtered-fungi.fa` and `data/pfam/PF00264.alignment.uniprot-cleaned-filtered-fungi-withoutgaps.fa`.

CD-HIT all kingdoms: `cd-hit -i data/pfam/PF00264.alignment.uniprot-cleaned-filtered-withoutgaps.fa -c 0.4 -n 2 -o data/pfam/PF00264.alignment.uniprot-cleaned-filtered-withoutgaps-cdhit0.4.fasta`

`cd-hit -i data/pfam/PF00264.alignment.uniprot-cleaned-filtered-withoutgaps.fa -c 0.5 -n 3 -o data/pfam/PF00264.alignment.uniprot-cleaned-filtered-withoutgaps-cdhit0.5.fasta`

CD-HIT fungi: `cd-hit -i data/pfam/PF00264.alignment.uniprot-cleaned-filtered-fungi-withoutgaps.fa -c 0.5 -n 3 -o data/pfam/PF00264.alignment.uniprot-cleaned-filtered-fungi-withoutgaps-cdhit0.5.fasta`

`cd-hit -i data/pfam/PF00264.alignment.uniprot-cleaned-filtered-fungi-withoutgaps.fa -c 0.6 -n 4 -o data/pfam/PF00264.alignment.uniprot-cleaned-filtered-fungi-withoutgaps-cdhit0.6.fasta`

To make the redundancy reduced alignment, run `python src/data-collection/make-redundancy-reduced-fasta.py`. This produces the files `data/pfam/PF00264.alignment.uniprot-cleaned-filtered-cdhit0.4.fasta` and `data/pfam/PF00264.alignment.uniprot-cleaned-filtered-cdhit0.5.fasta` and `data/pfam/PF00264.alignment.uniprot-cleaned-filtered-fungi-cdhit0.4.fasta` and `data/pfam/PF00264.alignment.uniprot-cleaned-filtered-fungi-cdhit0.5.fasta`.

Make trees for all kingdoms: `raxml-ng --msa data/pfam/PF00264.alignment.uniprot-cleaned-filtered-cdhit0.4.fasta --model JTT+G4 --prefix data/pfam/raxml/T8 --threads 7 --seed 2 --blopt nr_safe`

`raxml-ng --msa data/pfam/PF00264.alignment.uniprot-cleaned-filtered-cdhit0.5.fasta --model JTT+G4 --prefix data/pfam/raxml/T9 --threads 7 --seed 2 --blopt nr_safe`

Make the trees for fungi: `raxml-ng --msa data/pfam/PF00264.alignment.uniprot-cleaned-filtered-fungi-cdhit0.5.fasta --model JTT+G4 --prefix data/pfam/raxml/T10 --threads 7 --seed 2 --blopt nr_safe`

`raxml-ng --msa data/pfam/PF00264.alignment.uniprot-cleaned-filtered-fungi-cdhit0.6.fasta --model JTT+G4 --prefix data/pfam/raxml/T11 --threads 7 --seed 2 --blopt nr_safe`

To clean the trees, use regex "\.\d\/\d+-\d+".

Make trees for all kingdoms, cutoff 25: `raxml-ng --msa data/pfam/PF00264.alignment.uniprot-cleaned-filtered25-cdhit0.4.fasta --model JTT+G4 --prefix data/pfam/raxml/T12 --threads 7 --seed 2 --blopt nr_safe`

`raxml-ng --msa data/pfam/PF00264.alignment.uniprot-cleaned-filtered25-cdhit0.5.fasta --model JTT+G4 --prefix data/pfam/raxml/T13 --threads 7 --seed 2 --blopt nr_safe`

Make the trees for fungi, cutoff 25: `raxml-ng --msa data/pfam/PF00264.alignment.uniprot-cleaned-filtered25-fungi-cdhit0.5.fasta --model JTT+G4 --prefix data/pfam/raxml/T14 --threads 7 --seed 2 --blopt nr_safe`

`raxml-ng --msa data/pfam/PF00264.alignment.uniprot-cleaned-filtered25-fungi-cdhit0.6.fasta --model JTT+G4 --prefix data/pfam/raxml/T15 --threads 7 --seed 2 --blopt nr_safe`

To clean the trees, use regex "\.\d\/\d+-\d+".

#### Make itol files with cd-hit cluster size
To make itol label files with size of cd-hit clusters, run `src/itol-label-files/make-cdhit-cluster-size-label-files.py`.

#### Find representative sequences
HMMalign på seeds: `hmmalign --trim data/pfam/PF00264.hmm data/seeds-names.fa > data/seeds-names.hmmalign`.

Convert alignment to fasta and remove columns with dots: `seqconverter --remhmminsertcols -I stockholm -O fasta data/seeds-names.hmmalign > data/seeds-names.hmmalign.fa`

Remove gaps: `data/seeds-names.hmmalign-withoutgaps.fa`.

Add above sequences to the fastas for the tree: `data/pfam/PF00264.alignment.uniprot-cleaned-filtered-fungi-withoutgaps-withseeds.fa` and `data/pfam/PF00264.alignment.uniprot-cleaned-filtered-withoutgaps-withseeds.fa`

CD-HIT all kingdoms: `cd-hit -i data/pfam/PF00264.alignment.uniprot-cleaned-filtered-withoutgaps-withseeds.fa -c 0.4 -n 2 -o data/pfam/PF00264.alignment.uniprot-cleaned-filtered-withoutgaps-withseeds-cdhit0.4.fasta`

`cd-hit -i data/pfam/PF00264.alignment.uniprot-cleaned-filtered-withoutgaps-withseeds.fa -c 0.5 -n 3 -o data/pfam/PF00264.alignment.uniprot-cleaned-filtered-withoutgaps-withseeds-cdhit0.5.fasta`

CD-HIT fungi: `cd-hit -i data/pfam/PF00264.alignment.uniprot-cleaned-filtered-fungi-withoutgaps-withseeds.fa -c 0.5 -n 3 -o data/pfam/PF00264.alignment.uniprot-cleaned-filtered-fungi-withoutgaps-withseeds-cdhit0.5.fasta`

`cd-hit -i data/pfam/PF00264.alignment.uniprot-cleaned-filtered-fungi-withoutgaps-withseeds.fa -c 0.6 -n 4 -o data/pfam/PF00264.alignment.uniprot-cleaned-filtered-fungi-withoutgaps-withseeds-cdhit0.6.fasta`

Make itol files: `python src/itol-label-files/make-representatives-itol-label-files.py`.

##### For e-25 filtering
Add seeds to the fastas for the tree: `data/pfam/PF00264.alignment.uniprot-cleaned-filtered25-fungi-withoutgaps-withseeds.fa` and `data/pfam/PF00264.alignment.uniprot-cleaned-filtered25-withoutgaps-withseeds.fa`

CD-HIT all kingdoms: `cd-hit -i data/pfam/PF00264.alignment.uniprot-cleaned-filtered25-withoutgaps-withseeds.fa -c 0.4 -n 2 -o data/pfam/PF00264.alignment.uniprot-cleaned-filtered25-withoutgaps-withseeds-cdhit0.4.fasta`

`cd-hit -i data/pfam/PF00264.alignment.uniprot-cleaned-filtered25-withoutgaps-withseeds.fa -c 0.5 -n 3 -o data/pfam/PF00264.alignment.uniprot-cleaned-filtered25-withoutgaps-withseeds-cdhit0.5.fasta`

CD-HIT fungi: `cd-hit -i data/pfam/PF00264.alignment.uniprot-cleaned-filtered25-fungi-withoutgaps-withseeds.fa -c 0.5 -n 3 -o data/pfam/PF00264.alignment.uniprot-cleaned-filtered25-fungi-withoutgaps-withseeds-cdhit0.5.fasta`

`cd-hit -i data/pfam/PF00264.alignment.uniprot-cleaned-filtered25-fungi-withoutgaps-withseeds.fa -c 0.6 -n 4 -o data/pfam/PF00264.alignment.uniprot-cleaned-filtered25-fungi-withoutgaps-withseeds-cdhit0.6.fasta`

Make itol files: `python src/itol-label-files/make-representatives-itol-label-files.py`.

# Only fungi
Make fasta with only fungi: `python src/data-collection/make-fungi-fasta.py`. This creates the file `data/pfam/protein-matching-PF00264-fungi-shortheaders.fasta`.

`cd-hit -i data/pfam/protein-matching-PF00264-fungi-shortheaders.fasta -c 0.4 -n 2 -o data/pfam/protein-matching-PF00264-fungi-shortheaders-cdhit0.4.fasta`

HMMalign med trimming:
`hmmalign --trim data/pfam/PF00264.hmm data/pfam/protein-matching-PF00264-fungi-shortheaders-cdhit0.4-andseeds.fasta > data/pfam/protein-matching-PF00264-fungi-shortheaders-cdhit0.4-andseeds.fasta.hmmalign`

Convert to fasta online and make upper case: `protein-matching-PF00264-fungi-shortheaders-cdhit0.4-andseeds.fasta.hmmalign.fasta`

`raxml-ng --search1 --msa data/pfam/protein-matching-PF00264-fungi-shortheaders-cdhit0.4-andseeds.fasta.hmmalign.fasta --model JTT+G4 --prefix data/pfam/raxml/T4 --threads 8 --seed 2`

If it fails, run:
`raxml-ng --search1 --msa data/pfam/protein-matching-PF00264-fungi-shortheaders-cdhit0.4-andseeds.fasta.hmmalign.fasta --model JTT+G4 --prefix data/pfam/raxml/T4 --threads 8 --seed 2 --blopt nr_safe`

# Retrieve taxonomy
To make a tsv file with taxonomy for each uniprot hit, run: `python src/data-collection/make-taxonomy-file.py`. This creates the file `data/pfam/protein-matching-PF00264-fungi.tsv`.

# Find out which sequences represent the seeds
The fungal seeds are added to the reduced file: `data/pfam/protein-matching-PF00264-fungi-shortheaders-andseeds.fasta`.

Run CD-HIT: `cd-hit -i data/pfam/protein-matching-PF00264-fungi-shortheaders-andseeds.fasta -c 0.4 -n 2 -o data/pfam/protein-matching-PF00264-fungi-shortheaders-andseeds-cdhit0.4.fasta`

A file with which sequences represent each seed is manually created: `data/pfam/representative-sequences`.

An iTOL label file with this information is created manually: `data/itol-label-files/seed-representatives-not-filtered.txt`.

For the filtered ones, the equivalent precedure is used, and the final file is`data/itol-label-files/seed-representatives-filtered.txt`.

# Interproscan on fungal uniprot hits
(Delete this?) Transfer the hits file to the HPC: `scp data/pfam/protein-matching-PF00264-fungi-shortheaders-cdhit0.4.fasta idamei@transfer.gbar.dtu.dk:/work3/idamei/polyphenol-oxidases/interproscan-uniprot`

Copy the fasta to the HPC: `scp data/pfam/protein-matching-PF00264-fungi-shortheaders-cdhit0.4.fasta idamei@transfer.gbar.dtu.dk:/work3/idamei/polyphenol-oxidases/interproscan-uniprot/`

On the HPC:
Chunk fasta: `chunkfasta -c 20 -d polyphenol-oxidases/interproscan-uniprot polyphenol-oxidases/interproscan-uniprot/protein-matching-PF00264-fungi-shortheaders-cdhit0.4.fasta`

Run this:
`qrsh`
`/work3/idamei/bin/my_interproscan/interproscan-5.64-96.0/interproscan.sh -appl Pfam,SignalP_EUK,SignalP_GRAM_NEGATIVE,SignalP_GRAM_POSITIVE,Phobius -i /work3/idamei/polyphenol-oxidases/interproscan-uniprot/chunk00.fa -f tsv -o /work3/idamei/polyphenol-oxidases/interproscan-uniprot/chunk00.interproscan`

`/work3/idamei/bin/my_interproscan/interproscan-5.64-96.0/interproscan.sh -appl Pfam,SignalP_EUK,SignalP_GRAM_NEGATIVE,SignalP_GRAM_POSITIVE,Phobius -i /work3/idamei/polyphenol-oxidases/interproscan-uniprot/chunk01.fa -f tsv -o /work3/idamei/polyphenol-oxidases/interproscan-uniprot/chunk01.interproscan`
etc.

Download the results: `scp -r idamei@transfer.gbar.dtu.dk:/work3/idamei/polyphenol-oxidases/interproscan-uniprot data`

# Interproscan on uniprot hits for all kingdoms
Copy the fasta to the HPC: `scp data/pfam/protein-matching-PF00264-shortheaders-filtered-cdhit0.4.fasta idamei@transfer.gbar.dtu.dk:/work3/idamei/polyphenol-oxidases/interproscan-uniprot-allkingdoms`

On the HPC:
Chunk fasta: `chunkfasta -c 20 -d polyphenol-oxidases/interproscan-uniprot-allkingdoms polyphenol-oxidases/interproscan-uniprot-allkingdoms/protein-matching-PF00264-shortheaders-filtered-cdhit0.4.fasta`

Run this:
`qrsh`
`/work3/idamei/bin/my_interproscan/interproscan-5.64-96.0/interproscan.sh -appl Pfam,SignalP_EUK,SignalP_GRAM_NEGATIVE,SignalP_GRAM_POSITIVE -i /work3/idamei/polyphenol-oxidases/interproscan-uniprot-allkingdoms/chunk00.fa -f tsv -o /work3/idamei/polyphenol-oxidases/interproscan-uniprot-allkingdoms/chunk00.interproscan`

`/work3/idamei/bin/my_interproscan/interproscan-5.64-96.0/interproscan.sh -appl Pfam,SignalP_EUK,SignalP_GRAM_NEGATIVE,SignalP_GRAM_POSITIVE -i /work3/idamei/polyphenol-oxidases/interproscan-uniprot-allkingdoms/chunk01.fa -f tsv -o /work3/idamei/polyphenol-oxidases/interproscan-uniprot-allkingdoms/chunk01.interproscan`
etc.

Download the results: `scp -r idamei@transfer.gbar.dtu.dk:/work3/idamei/polyphenol-oxidases/interproscan-uniprot-allkingdoms data`

## Interproscan on proteome selected sequences
### All kingdoms, one per order
To make fasta with full-length sequences of the selected sequences from proteomes, run `python src/proteome-tree/get-full-length-selected-sequences.py`. This creates the file `data/proteome-tree/all-one_proteome_per_order_class.full-length.fa`.

Copy the fasta to the HPC: `scp data/proteome-tree/all-one_proteome_per_order.full-length.fa idamei@transfer.gbar.dtu.dk:/work3/idamei/polyphenol-oxidases/interproscan/proteomes-all-per-order-interproscan`

On the HPC:
Chunk fasta: `chunkfasta -c 20 -d polyphenol-oxidases/interproscan/proteomes-all-per-order-interproscan polyphenol-oxidases/interproscan/proteomes-all-per-order-interproscan/all-one_proteome_per_order.full-length.fa`

Run this:
`qrsh`
`/work3/idamei/bin/my_interproscan/interproscan-5.64-96.0/interproscan.sh -appl Pfam,SignalP_EUK,SignalP_GRAM_NEGATIVE,SignalP_GRAM_POSITIVE,Phobius -i /work3/idamei/polyphenol-oxidases/interproscan/proteomes-all-per-order-interproscan/chunk00.fa -f tsv -o /work3/idamei/polyphenol-oxidases/interproscan/proteomes-all-per-order-interproscan/chunk00.interproscan`

`/work3/idamei/bin/my_interproscan/interproscan-5.64-96.0/interproscan.sh -appl Pfam,SignalP_EUK,SignalP_GRAM_NEGATIVE,SignalP_GRAM_POSITIVE -i /work3/idamei/polyphenol-oxidases/interproscan/proteomes-all-per-order-interproscan/chunk01.fa -f tsv -o /work3/idamei/polyphenol-oxidases/interproscan/proteomes-all-per-order-interproscan/chunk01.interproscan`

etc.

Download the results: `scp -r idamei@transfer.gbar.dtu.dk:/work3/idamei/polyphenol-oxidases/interproscan/proteomes-all-per-order-interproscan data`

### Download domain data from interproscan
A json file was downloaded with all proteome data (e-mail correspondence uniprot): `data/proteome-tree/export.json`.

To download the domain data through the interproscan api, run `python src/proteome-tree/download-domain-data-api.py`. This downloads the files to `data/pfam/api`.

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
Files with proteins in each clade (determined from iTOL) were saved in `data/proteome-tree/clades`.

A csv file is made by running `python src/proteome-tree/clade-annotations.py`.

iTOL heatmap file is made by running `python src/itol-label-files/make-itol-label-files.py` and is shown on the species tree.

# MrBayes
On hal, go to the folder and start tmux: `tmux new -s mrbayes`

Start the job: `mpirun -np 6 mb gpu_run.nexus > log.txt`

Detach: control+b followed by d

To continue a mrbayes run, copy everything from the .ckp~ file into the bottom of the aignment nexus file. In the run nexus file, change the number of generations to the total number of generations you want and add append=yes, fx: `mcmc append=yes ngen=60000000 samplefreq=1000 nchains=3 file=out.nex`.

The run with all kingdoms, one per order was saved in `data/mrbayes/all`.

# Tmux
List sessions: `tmux list-sessions`
`tmux kill-session -t mysession`

# Fingerprint
The ids in each clade were saved in `data/mrbayes/all/clades`.

Fasta files are generated by running `python src/finger-print/make-clade-fastas.py`.