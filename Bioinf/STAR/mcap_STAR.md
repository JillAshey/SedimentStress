## STAR M. capitata

Date: 20220217

For sediment stress paper, we initially did not have Mcap metadata and so did not analyze Mcap samples. Now metadata is found so I need to re-run STAR using the Mcap version2 genomic info 

Deleted old Mcap genome directory off of server bc those files were all of old Mcap version. Now using files from http://cyanophora.rutgers.edu/montipora/ (version 2). Download and unzip all files before moving forward.

File quality and trimming has already been done for these samples, so will go right into STAR

a) After downloading onto server, unzip files 

```
ls
Montipora_capitata_HIv2.assembly.fasta.gz  Montipora_capitata_HIv2.genes.cds.fna.gz  Montipora_capitata_HIv2.genes.gff3.gz  Montipora_capitata_HIv2.genes.pep.faa.gz

gunzip *
```

Count number of seqs per file 

```
# Protein
zgrep -c ">" Montipora_capitata_HIv2.genes.pep.faa 
55086

# Transcript
zgrep -c ">" Montipora_capitata_HIv2.genes.cds.fna 
55086
```

b) Generate genome index

STAR needs transcript_id identifier to run properly, so added that identifier for Pacuta in R. Renamed gff file as Mcap.gff.annotations.fixed_transcript.gff3 and uploaded to server.

Moved all old STAR output into new directory called 


```
cd /data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR
mkdir Mcap_V1
mv *mcap* Mcap_V1/
```

Run genome index step of STAR

Write as a script and submit as job 

```
mkdir GenomeIndex_mcapV2

nano star_index_mcap.sh

#!/bin/bash
#SBATCH -t 48:00:00
#SBATCH --nodes=1 --ntasks-per-node=2
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Francois_data/Hawaii/scripts
#SBATCH --error="star_index_mcap_out_error"
#SBATCH --output="star_index_mcap_out"

module load STAR/2.5.3a-foss-2016b # on bluewaves

STAR --runThreadN 10 --runMode genomeGenerate --genomeDir /data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/GenomeIndex_mcapV2 --genomeFastaFiles /data/putnamlab/jillashey/genome/Mcap/Montipora_capitata_HIv2.assembly.fasta --sjdbGTFfile /data/putnamlab/jillashey/genome/Mcap/Mcap.gff.annotations.fixed_transcript.gff3

sbatch star_index_mcap.sh
```

Submitted batch job 1978534
