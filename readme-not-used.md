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

# Filter uniprot members
Filter fasta `python src/data-collection/filter.py`. This generates the file `data/pfam/protein-matching-PF00264-shortheaders-filtered25.fasta`.