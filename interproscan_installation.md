# Installation of Interproscan on the HPC
## Installing Interproscan (see https://interproscan-docs.readthedocs.io/en/latest/Introduction.html)
### Downloading interproscan
`mkdir /work3/idamei/bin/my_interproscan/`
`cd /work3/idamei/bin/my_interproscan`
`wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.64-96.0/interproscan-5.64-96.0-64-bit.tar.gz`
`wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.64-96.0/interproscan-5.64-96.0-64-bit.tar.gzmd`

### checksum to confirm the download was successful:
`md5sum -c interproscan-5.6-96.0-64-bit.tar.gz.md`

Must return *interproscan-5.64-96.0-64-bit.tar.gz: OK*

### Extract the tar ball:
`tar -pxvzf interproscan-5.64-96.0-*-bit.tar.gz`

### Index hmm models
`python3 setup.py -f interproscan.properties`

## Installing Java
The JDK build was downloaded at https://jdk.java.net/21/ (the Linux/x64 version was chosen) and transferred to `/work3/idamei/bin`.

And unpacked:
`cd /work3/idamei/bin`
`tar xvf openjdk-13*_bin.tar.gz`

## Installing pftools (see https://github.com/sib-swiss/pftools3/blob/master/INSTALL)
`cd /work3/idamei/bin`

`git clone https://github.com/sib-swiss/pftools3.git`
`cd pftools3`
`mkdir build`
`cd build/`
`cmake -DUSE_AFFINITY=OFF -DUSE_PCRE2=OFF -DUSE_PCRE=ON -DCMAKE_INSTALL_PREFIX:PATH=/work3/idamei/bin/pft ..`
`make`
`make install`
`make test`

## Adding pftools path to interproscan.properties
In `/work3/idamei/bin/my_interproscan/interproscan-5.64-96.0/interproscan.properties`, the following two files were added:

`binary.prosite.pfscanv3.path=/work3/idamei/bin/pft/bin/pfscanV3`
`binary.prosite.pfscanv3.path=/work3/idamei/bin/pft/bin/pfsearchV3`

## Interproscan is run this way:
`/work3/idamei/bin/my_interproscan/interproscan-5.64-96.0/interproscan.sh -i [input_fasta_filename] -f tsv -o [output_filename]`