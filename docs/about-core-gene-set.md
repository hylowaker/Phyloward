### (Document Work in Progress)

## Establishment of archaeal core gene set

The 74 archaeal core gene set was determined by the following method. 
Complete archaeal genome sequences available from the [EzBioCloud](https://www.ezbiocloud.net) database were used for further steps. 
Among all archaeal genomes in the database, single complete genome with best quality was selected for each species to avoid contaminated genomes and the bias of phylum among the member genomes. 
Total 225 genomes were selected as representatives, and their protein coding genes were used as the target by hmmsearch software in HMMER package. 
With every hidden Markov model (HMM) profile in [TIGRFAMs](http://tigrfams.jcvi.org/cgi-bin/index.cgi) database and [Pfam-A](https://pfam.xfam.org) database, hmmsearch was performed on all protein coding genes of the 225 genomes. 
The gene presented as single-copy in at least 95% of 225 genomes is considered as core gene. Total 71 genes were selected as core gene on this step.