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

# CD-HIT, raxml (not used)
Filter fasta `python src/data-collection/filter.py`. This generates the file `data/pfam/protein-matching-PF00264-shortheaders-filtered.fasta`.

CD-HIT: `cd-hit -i data/pfam/protein-matching-PF00264-shortheaders-filtered.fasta -c 0.4 -n 2 -o data/pfam/protein-matching-PF00264-shortheaders-filtered-cdhit0.4.fasta`

HMMalign: `hmmalign --trim data/pfam/PF00264.hmm data/pfam/protein-matching-PF00264-shortheaders-filtered-cdhit0.4.fasta > data/pfam/protein-matching-PF00264-shortheaders-filtered-cdhit0.4.fasta.hmmaligntrim`.

Convert output file to fasta online and make upper case (http://sequenceconversion.bugaco.com/converter/biology/sequences/stockholm_to_fasta.php): `data/pfam/protein-matching-PF00264-shortheaders-filtered-cdhit0.4.fasta.hmmaligntrim.fasta`

Make tree: `raxml-ng --search1 --msa data/pfam/protein-matching-PF00264-shortheaders-filtered-cdhit0.4.fasta.hmmaligntrim.fasta --model JTT+G4 --prefix data/pfam/raxml/T5 --threads 8 --seed 2`

# Making trimmed fasta from alignment (not used)
The alignment is converted to fasta: `seqconverter -I stockholm -O fasta data/pfam/PF00264.alignment.uniprot > data/pfam/PF00264.alignment.uniprot.fa`. The output file is copied to `data/pfam/PF00264-trimmed.fa` where gaps are removed. 

The alignment is cleaned from columns with insertions by running python `src/data-collection/clean-hmm-alignment.py`. This created the fasta file `data/pfam/PF00264.alignment.uniprot-cleaned.fa`.

A trimmed fasta was made: `data/pfam/PF00264.alignment.uniprot-nogaps.fa` by removing all gaps from the alignment.

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

# Make itol files with cd-hit cluster size
To make itol label files with size of cd-hit clusters, run `src/itol-label-files/make-cdhit-cluster-size-label-files.py`.

#### Find representative sequences
HMMalign på seeds: `hmmalign --trim data/pfam/PF00264.hmm data/seeds-names.fa > data/seeds-names.hmmalign`.

Convert alignment to fasta and remove columns with dots: `seqconverter --remhmminsertcols -I stockholm -O fasta data/seeds-names.hmmalign > data/seeds-names.hmmalign.fa`

(NEW) `seqconverter -I stockholm -O fasta data/seeds-names.hmmalign > data/seeds-names.hmmalign.fa`. Copy this file to `data/seeds-names-trimmed.fa` and remove gaps.

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

# Find out which sequences represent the seeds
The fungal seeds are added to the reduced file: `data/pfam/protein-matching-PF00264-fungi-shortheaders-andseeds.fasta`.

Run CD-HIT: `cd-hit -i data/pfam/protein-matching-PF00264-fungi-shortheaders-andseeds.fasta -c 0.4 -n 2 -o data/pfam/protein-matching-PF00264-fungi-shortheaders-andseeds-cdhit0.4.fasta`

A file with which sequences represent each seed is manually created: `data/pfam/representative-sequences`.

An iTOL label file with this information is created manually: `data/itol-label-files/seed-representatives-not-filtered.txt`.

For the filtered ones, the equivalent precedure is used, and the final file is`data/itol-label-files/seed-representatives-filtered.txt`.

# Only fungi
Make fasta with only fungi: `python src/data-collection/make-fungi-fasta.py`. This creates the file `data/pfam/protein-matching-PF00264-fungi-shortheaders.fasta`.

`cd-hit -i data/pfam/protein-matching-PF00264-fungi-shortheaders.fasta -c 0.4 -n 2 -o data/pfam/protein-matching-PF00264-fungi-shortheaders-cdhit0.4.fasta`

HMMalign med trimming:
`hmmalign --trim data/pfam/PF00264.hmm data/pfam/protein-matching-PF00264-fungi-shortheaders-cdhit0.4-andseeds.fasta > data/pfam/protein-matching-PF00264-fungi-shortheaders-cdhit0.4-andseeds.fasta.hmmalign`

Convert to fasta online and make upper case: `protein-matching-PF00264-fungi-shortheaders-cdhit0.4-andseeds.fasta.hmmalign.fasta`

`raxml-ng --search1 --msa data/pfam/protein-matching-PF00264-fungi-shortheaders-cdhit0.4-andseeds.fasta.hmmalign.fasta --model JTT+G4 --prefix data/pfam/raxml/T4 --threads 8 --seed 2 --blopt nr_safe`

# Aclust tree
Aclust tree with full length was made on the HPC and saved in `data/trees/aclust-uniprot-al-kingdoms-filtered.nwk`. This file was used as input: `data/pfam/protein-matching-PF00264-shortheaders-filtered-cdhit0.4.fasta`.

An aclust tree with trimmed sequences was made on the HPC and saved in `data/trees/aclust-uniprot-al-kingdoms-filtered.nwk`. This file was used as input: `data/pfam/protein-matching-PF00264-shortheaders-filtered-cdhit0.4.fasta.hmmaligntrim.fasta.aclustinput`.

