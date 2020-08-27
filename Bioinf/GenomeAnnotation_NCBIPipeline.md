### Annotating genomes with NCBI pipeline


Pipeline significance

- Produces annotation of coding and non-coding genes, transcripts and proteins for finished and unfinished genome assemblies
- Provides content for other NCBI resources (Nucleotide, Protein, BLAST, Gene, Map Viewer)

![](https://www.ncbi.nlm.nih.gov/core/assets/genome/images/Pipeline_RFAM.png)


#### 1) Retrieving inputs

Inputs will be: transcript file, protein file, genome file

#### 2) Genome seq. masking

Using RepeatMasker

#### 3) Alignment

Transcript alignment uses BLAST and Splign

###### Need to look into RefSeq more to see if P. dam genomic info is stored there. May change some settings for analysis

Protein alignment uses BLAST and ProSplign

#### 4) Gene prediction

Alignments are passed to Gnomon for gene prediction. Gnomon is a two-step process

- Chainer - assembles partial alignments into putative gene models 
- *ab initio* predictions

First predictions of Gnomon is refined by alignment to protein seq. database. Those alignments are added, and Gnomon is run again. 

#### 5) Choosing best model

Final model will be based on Gnomon annotations (and possibly RefSeq sequences?). Annotations using RefSeq and Gnomon can be integrated. Gnomon predictions will be labelled with prefixes (XM_ protein coding transcripts, XR_ non-coding transcripts, XP_ proteins) to distinguish results from RefSeq (NM_ protein coding transcripts, NR_ non-coding transcripts, NP_ proteins)







