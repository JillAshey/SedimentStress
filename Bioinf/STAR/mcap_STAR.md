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

20220219 - okay so it ran, but something got messed up. there is no info in the genome index for sjdbList files. It might be because of how I separated out the transcript_id in R? Going to go back now and check code

20220221 - I separated out the transcript_id in R using a different method, copied the new gff fixed file to server, and now going to re-run ```star_index_mcap.sh```

Submitted batch job 1979176

It seems to have worked! Good amount of info in the sjdbList files. Now to align reads against mcap

c) Align reads 

The mcap samples are: `5_2.fastq.gz, 10_2.fastq.gz, 13_2.fastq.gz, 14_2.fastq.gz, 15_2.fastq.gz, 16_2.fastq.gz, 19_2.fastq.gz, 20_2.fastq.gz, 40_2.fastq.gz, 45_2.fastq.gz, 48_2.fastq.gz`. Going to separate them out to run on STAR first and then will run all other samples


Link to trimmed files 

```
mkdir AlignReads_mcapV2
cd AlignReads_mcapV2
ln -s /data/putnamlab/jillashey/Francois_data/Hawaii/data/trimmed/*trim.fq .
```

Move mcap samples to their own directory 

```
mkdir mcap_only
mv 5_2.fastq.trim.fq 10_2.fastq.trim.fq 13_2.fastq.trim.fq 14_2.fastq.trim.fq 15_2.fastq.trim.fq 16_2.fastq.trim.fq 19_2.fastq.trim.fq 20_2.fastq.trim.fq 40_2.fastq.trim.fq 45_2.fastq.trim.fq 48_2.fastq.trim.fq mcap_only/
cd mcap_only
```

Let's do a test sample first 

```
nano test_5_2_align.sh

#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="test_align_mcap_out_error"
#SBATCH --output="test_align_mcap_out_error"

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads --quantMode TranscriptomeSAM --outTmpDir /data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/test_TMP --readFilesIn 5_2.fastq.trim.fq --genomeDir /data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/GenomeIndex_mcapV2/ --twopassMode Basic --twopass1readsN -1 --outStd Log BAM_Unsorted BAM_Quant --outSAMtype BAM Unsorted SortedByCoordinate --outReadsUnmapped Fastx --outFileNamePrefix 5_2_test.

sbatch test_5_2_align.sh
```

Submitted batch job 1979184

Keep getting this error: EXITING because of fatal ERROR: could not make temporary directory: /data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/
SOLUTION: (i) please check the path and writing permissions 
 (ii) if you specified --outTmpDir, and this directory exists - please remove it before running STAR

Feb 21 23:43:57 ...... FATAL ERROR, exiting

Not sure why... going to try later 

Okay tried again. I had the --outTmpDir going to `/data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/AlignReads_mcapV2` but changed it so it would go to `/data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/test_TMP `. Submitted batch job 1979524. Success! Took ~75 mins 

Now will run w/ all mcap samples 

```
nano AlignReads_mcap_only.sh

#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="Align_mcap_only_out_error"
#SBATCH --output="Align_mcap_only_out"

module load STAR/2.5.3a-foss-2016b

F=/data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/AlignReads_mcapV2/mcap_only

array1=($(ls $F/*trim.fq))
for i in ${array1[@]}
do
STAR --runMode alignReads --quantMode TranscriptomeSAM --outTmpDir ${i}_TMP --readFilesIn ${i} --genomeDir /data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/GenomeIndex_mcapV2/ --twopassMode Basic --twopass1readsN -1 --outStd Log BAM_Unsorted BAM_Quant --outSAMtype BAM Unsorted SortedByCoordinate --outReadsUnmapped Fastx --outFileNamePrefix ${i}.
done 

sbatch --array=1-11 AlignReads_mcap_only.sh # submit as an array job 
```

Submitted batch job 1979526
