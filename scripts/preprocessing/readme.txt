## Preprocessing pipeline of circulating microRNA data ##

These scripts were used to process the Illumina sequenced small-RNAs (miRNAs). Preprocessing was conducted using HPC provided by CSC (Puhti, https://www.puhti.csc.fi/public/).

Scripts will perform:
  - Conversion of RNA into DNA
  - Indexing genome for alignment (Bowtie)
  - Preprocessing of the data
  - Merging data prior alignment
  - Alignment with Bowtie
  - Count matrix creation from aligned reads
