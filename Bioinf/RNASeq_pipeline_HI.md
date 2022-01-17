## Pipeline for HI species - sediment stress 

Species: _Pocillopora acuta_ and _Porites lobata_

Note: some of the paths may not be correct, as I have made new directories and moved files around since running this pipeline.

### 1) Check file integrity 

a) Count all files to make sure all downloaded

```
ls -1 | wc -l
```

b) Verify data integrity with md5sum

```
md5sum *.fastq > checkmd5.md5
md5sum -c checkmd5.md5
```
Should output 'OK' next to each file name 

c) Count number of reads per file 

Some files have different @___. There are: HISEQ, HWI

```
zgrep -c "HISEQ" *.fastq
zgrep -c "HWI" *.fastq
```

xxxxxxxxxx

### 2) Run FastQC

a) Make folders for raw FastQC results and scripts
```
cd Francois_data/Hawaii
mkdir fastqc_results/raw
mkdir scripts
cd scripts
```

b) Write script for checking quality with FastQC and submit as job on bluewaves

```
nano fastqc_raw.sh 

#!/bin/bash
#SBATCH -t 18:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@my.uri.edu
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Francois_data/Hawaii/scripts
#SBATCH --error="fastqc_out_raw_error"
#SBATCH --output="fastqc_out_raw"

module load FastQC/0.11.8-Java-1.8 #another version in Bluewaves if this one doesn't work

for file in /data/putnamlab/jillashey/Francois_data/Hawaii/data/raw/*.fastq.gz
do
fastqc $file --outdir /data/putnamlab/jillashey/Francois_data/Hawaii/fastqc_results/raw/
done

sbatch fastqc_raw.sh 
```

b) Run MultiQC. Pretty fast, so don't need to submit job for it 

```
module load MultiQC/1.7-foss-2018b-Python-2.7.15
multiqc /data/putnamlab/jillashey/Francois_data/Hawaii/fastqc_results/raw/*fastqc.zip -o /data/putnamlab/jillashey/Francois_data/Hawaii/multiqc_results/raw
```

c) Copy MultiQC files to local computer

Add multiQC plots here 

### 4) Trim reads with Trimmomatic 

a) Make trimmed reads folder in all other results folders 

```
cd Francois_data/Hawaii
mkdir data/trimmed fastqc_results/trimmed multiqc_results/trimmed
```

b) If necessary, unzip fastqc files, as Trimmomatic can't process zipped files

```
nano trimmomatic.sh

#!/bin/bash
#SBATCH -t 48:00:00
#SBATCH --nodes=1 --ntasks-per-node=2
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@my.uri.edu
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Francois_data/Hawaii/scripts
#SBATCH --error="trimmomatic_out_error"
#SBATCH --output="trimmomatic_out"

module load Trimmomatic/0.38-Java-1.8

for file in /data/putnamlab/jillashey/Francois_data/Hawaii/data/raw/*.fastq
do
        java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.38.jar SE -phred33 $file $file.trim.fq ILLUMINACLIP:/data/putnamlab/jillashey/Francois_data/Hawaii/data/Illumina_adapter_reads_PE_SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20 >> amtTrimmed.txt
done

sbatch trimmomatic.sh
```

d) Move trimmed files to their own folder 

```
mv *trim.fq /data/putnamlab/jillashey/Francois_data/Hawaii/data/trimmed
```

### 5) Check quality of trimmed files 

a) In trimmed read folder, check number of files 

```
ls -1 | wc -l
```

b) Check number of reads

```
zgrep -c "HISEQ" *trim.fq
zgrep -c "HWI" *trim.fq
```

c) Run FastQC on trimmed data

```
nano fastqc_trim.sh

#!/bin/bash
#SBATCH -t 48:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=150GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@my.uri.edu
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Francois_data/Hawaii/scripts
#SBATCH --error="fastqc_out_trim_error"
#SBATCH --output="fastqc_out_trim"

module load FastQC/0.11.8-Java-1.8 #another version in Bluewaves if this one doesn't work

for file in /data/putnamlab/jillashey/Francois_data/Hawaii/data/trimmed/*.trim.fq.gz
do
fastqc $file --outdir /data/putnamlab/jillashey/Francois_data/Hawaii/fastqc_results/trimmed/
done

sbatch fastqc_trim.sh
```

d) Run MultiQC on trimmed data and copy files to local computer

```
module load MultiQC/1.7-foss-2018b-Python-2.7.15
multiqc /data/putnamlab/jillashey/Francois_data/Hawaii/fastqc_results/trimmed/*fastqc.zip -o /data/putnamlab/jillashey/Francois_data/Hawaii/multiqc_results/trimmed
```

Add multiQC plots here

### 6) Align reads with STAR

Before running STAR, I added the identifier 'transcript_id=' to the last column of the gff file in R. STAR needs this identifier to run and most of the gffs I used don't have it. Code to edit gff files herexxxxxxx

#### P. acuta

a) Generate genome index

```
mkdir /data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/GenomeIndex_pacuta

module load STAR/2.5.3a-foss-2016b # on bluewaves

STAR --runThreadN 10 --runMode genomeGenerate --genomeDir /data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/GenomeIndex_pacuta --genomeFastaFiles /data/putnamlab/jillashey/genome/Pacuta/Pocillopora_acuta_HIv1.assembly.purged.fasta --sjdbGTFfile /data/putnamlab/jillashey/genome/Pacuta/Pacuta.gff.annotations.fixed_transcript.gff3
```

Said it was starting, but I stopped it after 20 mins because it shouldn't be taking that long. There was nothing in the GenomeIndex_pacuta directory either. 

```
Jan 10 13:30:46 ..... started STAR run
Jan 10 13:30:47 ... starting to generate Genome files
^C
```

Write as a script and submit as job 

```
nano star_index_pacuta.sh

#!/bin/bash
#SBATCH -t 48:00:00
#SBATCH --nodes=1 --ntasks-per-node=2
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Francois_data/Hawaii/scripts
#SBATCH --error="trimmomatic_out_error"
#SBATCH --output="trimmomatic_out"

module load STAR/2.5.3a-foss-2016b # on bluewaves

STAR --runThreadN 10 --runMode genomeGenerate --genomeDir /data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/GenomeIndex_pacuta --genomeFastaFiles /data/putnamlab/jillashey/genome/Pacuta/Pocillopora_acuta_HIv1.assembly.purged.fasta --sjdbGTFfile /data/putnamlab/jillashey/genome/Pacuta/Pacuta.gff.annotations.fixed_transcript.gff3

sbatch star_index_pacuta.sh
```

Submitted batch job 1960631. Worked! Took ~45 mins

b) Align reads to genome 

Link to trimmed files 

```
mkdir AlignReads_pacuta
cd AlignReads_pacuta
ln -s /data/putnamlab/jillashey/Francois_data/Hawaii/data/trimmed/*trim.fq .
```

Let's do a test sample first 

```
nano test_11_2_align.sh

#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="test_align_pacuta_out_error"
#SBATCH --output="test_align_pacuta_out_error"

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads --quantMode TranscriptomeSAM --outTmpDir /data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/test_TMP --readFilesIn 11_2.fastq.trim.fq --genomeDir /data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/GenomeIndex_pacuta/ --twopassMode Basic --twopass1readsN -1 --outStd Log BAM_Unsorted BAM_Quant --outSAMtype BAM Unsorted SortedByCoordinate --outReadsUnmapped Fastx --outFileNamePrefix 11_2_test.

sbatch test_11_2_align.sh
```

Submitted batch job 1960632

Moving all the Pacuta samples into their own directory and will run those first. Then I can go back and run the other species samples against the Pacuta genome.

```
mkdir pacuta_only

# move pacuta samples to new directory 
mv 1_2.fastq.trim.fq 2_2.fastq.trim.fq 4_2.fastq.trim.fq 28_2.fastq.trim.fq 35_2.fastq.trim.fq 36_2.fastq.trim.fq 38_2.fastq.trim.fq 39_2.fastq.trim.fq 41_2.fastq.trim.fq 42_2.fastq.trim.fq 47_2.fastq.trim.fq pacuta_only/

# did not move 11_2 yet because it is getting analyzed currently as the test run
```

Success! Took ~1.5 hrs

Now will run w/ all pacuta samples 

```
nano AlignReads_pacuta_only.sh

#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="Align_pacuta_only_out_error"
#SBATCH --output="Align_pacuta_only_out"

module load STAR/2.5.3a-foss-2016b

F=/data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/AlignReads_pacuta/pacuta_only

array1=($(ls $F/*trim.fq))
for i in ${array1[@]}
do
STAR --runMode alignReads --quantMode TranscriptomeSAM --outTmpDir ${i}_TMP --readFilesIn ${i} --genomeDir /data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/GenomeIndex_pacuta --twopassMode Basic --twopass1readsN -1 --outStd Log BAM_Unsorted BAM_Quant --outSAMtype BAM Unsorted SortedByCoordinate --outReadsUnmapped Fastx --outFileNamePrefix ${i}.
done 

sbatch --array=1-11 AlignReads_pacuta_only.sh # submit as an array job 
```

Submitted batch job 1960633

Told me job was done, took ~15 hours to run. for some reason, no output for 4_2 

Align 4_2 individually 

```
nano 4_2_align.sh

#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="4_2_align_pacuta_out_error"
#SBATCH --output="4_2_align_pacuta_out"

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads --quantMode TranscriptomeSAM --outTmpDir /data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/test_TMP --readFilesIn 4_2.fastq.trim.fq --genomeDir /data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/GenomeIndex_pacuta/ --twopassMode Basic --twopass1readsN -1 --outStd Log BAM_Unsorted BAM_Quant --outSAMtype BAM Unsorted SortedByCoordinate --outReadsUnmapped Fastx --outFileNamePrefix 4_2.

sbatch 4_2_align.sh
```

Submitted batch job 1960686. cancelled job because it was pending for >1 day and resubmitted -- success!

Will now align the rest of the samples against pacuta genome 

```
nano AlignReads_pacuta.sh

#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="Align_pacuta_out_error"
#SBATCH --output="Align_pacuta_out"

module load STAR/2.5.3a-foss-2016b

F=/data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/AlignReads_pacuta

array1=($(ls $F/*trim.fq))
for i in ${array1[@]}
do
STAR --runMode alignReads --quantMode TranscriptomeSAM --outTmpDir ${i}_TMP --readFilesIn ${i} --genomeDir /data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/GenomeIndex_pacuta --twopassMode Basic --twopass1readsN -1 --outStd Log BAM_Unsorted BAM_Quant --outSAMtype BAM Unsorted SortedByCoordinate --outReadsUnmapped Fastx --outFileNamePrefix ${i}.
done 

sbatch --array=1-52%15 AlignReads_pacuta.sh # submit as an array job 
```

Submitted batch job 1960658

#### P. lobata 

a) Generate genome index 

Using P.lutea genome to analyze p.lobata in this study

Attempt 1

```
module load STAR/2.5.3a-foss-2016b

STAR --runThreadN 10 --runMode genomeGenerate --genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_plob/ --genomeFastaFiles /data/putnamlab/jillashey/genome/Plutea/plut_final_2.1.fasta --sjdbGTFfile /data/putnamlab/jillashey/genome/Plutea/plut2v1.1.genes.gff3 
```
not working at the moment, may be something to do with fasta/gff file labels in the files themselves 
error:
```
Jul 30 16:41:17 ..... started STAR run
Jul 30 16:41:17 ... starting to generate Genome files
Jul 30 16:41:31 ... starting to sort Suffix Array. This may take a long time...
Jul 30 16:41:35 ... sorting Suffix Array chunks and saving them to disk...
Jul 30 16:43:54 ... loading chunks from disk, packing SA...
Jul 30 16:45:54 ... finished generating suffix array
Jul 30 16:45:54 ... generating Suffix Array index
Jul 30 16:47:54 ... completed Suffix Array index
Jul 30 16:47:54 ..... processing annotations GTF
terminate called after throwing an instance of 'std::out_of_range'
  what():  vector::_M_range_check: __n (which is 0) >= this->size() (which is 0)
Aborted (core dumped)
```
Copied the files from the tufts server to Desktop to Bluewaves to see if they will work

Attempt 2

```
module load STAR/2.5.3a-foss-2016b

STAR --runThreadN 10 --runMode genomeGenerate --genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_plob/ --genomeFastaFiles /data/putnamlab/jillashey/genome/plut_final_2.1.fasta --sjdbGTFfile /data/putnamlab/jillashey/genome/plut2v1.1.genes.gff3 
```

Attempt 3
STAR needs transcript ID to generate genome index. In R, I replaced ID with transcript_id in 9th row and generated an updated gff file. 

```
# getting new gff file from local computer
scp /Users/jillashey/Desktop/Plut.GFFannotation.fixed.gff jillashey@bluewaves.uri.edu:/data/putnamlab/jillashey/genome/Plutea/
```

```
module load STAR/2.5.3a-foss-2016b

STAR --runThreadN 10 --runMode genomeGenerate --genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_plob/ --genomeFastaFiles /data/putnamlab/jillashey/genome/plut_final_2.1.fasta --sjdbGTFfile /data/putnamlab/jillashey/genome/Plutea/Plut.GFFannotation.fixed.gff
```

This worked!!

Attempt 4
STAR needs transcript ID to generate genome index. In R, we removed some info from ID and added it to transcript id, which was added to gene column

```
# getting new gff file from local computer
scp /Users/jillashey/Desktop/Plut.GFFannotation.fixed_transcript.gff jillashey@bluewaves.uri.edu:/data/putnamlab/jillashey/genome/Plutea/
```

```
module load STAR/2.5.3a-foss-2016b

STAR --runThreadN 10 --runMode genomeGenerate --genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_plob/ --genomeFastaFiles /data/putnamlab/jillashey/genome/Plutea/plut_final_2.1.fasta --sjdbGTFfile /data/putnamlab/jillashey/genome/Plutea/Plut.GFFannotation.fixed_transcript.gff
```

Worked!

b) Align reads to genome 

#### P. lobata (using P. lutea genome)

a) make directories for STAR output
```
```


b) Obtain reference genome (.fasta or .fna) and reference annotation (.gtf or .gff). Unzip these files

```
# P.lutea genome info
/data/putnamlab/REFS/Plutea
# plut_final_2.1.fasta.gz -- genome file 
# plut2v1.1.genes.gff3.gz -- annotation file 
# Move to my own genome folder because REFS isn't allowing me to work with files 

# Plutea unzip
gunzip plut_final_2.1.fasta.gz
gunzip plut2v1.1.genes.gff3.gz
```

c) Unzip fastq files so they can properly analyzed by STAR

```
#!/bin/bash
#SBATCH -t 12:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@my.uri.edu
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Francois_data/scripts
#SBATCH --error="gunzip_trim_out_error"
#SBATCH --output="gunzip_trim_out"

for file in /data/putnamlab/jillashey/Francois_data/data/trimmed/*.fq.gz
do
gunzip $file 
done
```

d) Generate genome index 
For p.lobata (using P.lutea genome)

Attempt 1

```
module load STAR/2.5.3a-foss-2016b

STAR --runThreadN 10 --runMode genomeGenerate --genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_plob/ --genomeFastaFiles /data/putnamlab/jillashey/genome/Plutea/plut_final_2.1.fasta --sjdbGTFfile /data/putnamlab/jillashey/genome/Plutea/plut2v1.1.genes.gff3 
```
not working at the moment, may be something to do with fasta/gff file labels in the files themselves 
error:
```
Jul 30 16:41:17 ..... started STAR run
Jul 30 16:41:17 ... starting to generate Genome files
Jul 30 16:41:31 ... starting to sort Suffix Array. This may take a long time...
Jul 30 16:41:35 ... sorting Suffix Array chunks and saving them to disk...
Jul 30 16:43:54 ... loading chunks from disk, packing SA...
Jul 30 16:45:54 ... finished generating suffix array
Jul 30 16:45:54 ... generating Suffix Array index
Jul 30 16:47:54 ... completed Suffix Array index
Jul 30 16:47:54 ..... processing annotations GTF
terminate called after throwing an instance of 'std::out_of_range'
  what():  vector::_M_range_check: __n (which is 0) >= this->size() (which is 0)
Aborted (core dumped)
```
Copied the files from the tufts server to Desktop to Bluewaves to see if they will work

Attempt 2

```
module load STAR/2.5.3a-foss-2016b

STAR --runThreadN 10 --runMode genomeGenerate --genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_plob/ --genomeFastaFiles /data/putnamlab/jillashey/genome/plut_final_2.1.fasta --sjdbGTFfile /data/putnamlab/jillashey/genome/plut2v1.1.genes.gff3 
```

Attempt 3
STAR needs transcript ID to generate genome index. In R, I replaced ID with transcript_id in 9th row and generated an updated gff file. 

```
# getting new gff file from local computer
scp /Users/jillashey/Desktop/Plut.GFFannotation.fixed.gff jillashey@bluewaves.uri.edu:/data/putnamlab/jillashey/genome/Plutea/
```

```
module load STAR/2.5.3a-foss-2016b

STAR --runThreadN 10 --runMode genomeGenerate --genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_plob/ --genomeFastaFiles /data/putnamlab/jillashey/genome/plut_final_2.1.fasta --sjdbGTFfile /data/putnamlab/jillashey/genome/Plutea/Plut.GFFannotation.fixed.gff
```

This worked!!

Attempt 4
STAR needs transcript ID to generate genome index. In R, we removed some info from ID and added it to transcript id, which was added to gene column

```
# getting new gff file from local computer
scp /Users/jillashey/Desktop/Plut.GFFannotation.fixed_transcript.gff jillashey@bluewaves.uri.edu:/data/putnamlab/jillashey/genome/Plutea/
```

```
module load STAR/2.5.3a-foss-2016b

STAR --runThreadN 10 --runMode genomeGenerate --genomeDir /data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/GenomeIndex_plob/ --genomeFastaFiles /data/putnamlab/jillashey/genome/Plutea/plut_final_2.1.fasta --sjdbGTFfile /data/putnamlab/jillashey/genome/Plutea/Plut.GFFannotation.fixed_transcript.gff
```

Worked!


e) Align reads

For P.lob (using p.lutea genome)

Try with one sample first 

```
nano star_align_plob_test.sh
```

```
#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/test_6_1_TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/6_1.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_plob/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Plutea/Plut.GFFannotation.fixed_transcript.gff
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_plob/test_6_1.
```

```
sbatch star_align_plob_test.sh
```

Submitted batch job 1667173

Won't working...no log.out.final

Attempt 1)

With all 64 samples 

```
nano star_align_batch_script_plob.sh
```

```
#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

FILENAME=$(ls -l /data/putnamlab/jillashey/Francois_data/data/trimmed/ | sed -n $((${SLURM_ARRAY_TASK_ID}+1))p | awk '{print $9}' )

STAR --runMode alignReads --quantMode TranscriptomeSAM --outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/${FILENAME}_TMP --readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/${FILENAME} --genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_plob/ --sjdbGTFfeatureExon exon --sjdbGTFfile /data/putnamlab/jillashey/genome/Plutea/Plut.GFFannotation.fixed.gff --sjdbGTFtagExonParentTranscript ID --twopassMode Basic --twopass1readsN -1 --outStd Log BAM_Unsorted BAM_Quant --outSAMtype BAM Unsorted SortedByCoordinate --outReadsUnmapped Fastx --outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_plob/${FILENAME}.
```

```
# Submit as an array job 
sbatch --array 1-64 star_align_batch_script_plob.sh
```
Submitted batch job 1667295

in Log.out, there is this warning: WARNING: while processing sjdbGTFfile=/data/putnamlab/jillashey/genome/Plutea/Plut.GFFannotation.fixed.gff: no gene_id for line: - all lines seem to be exons

No log.final.out, meaning it did not finish the complete alignment. 

Looking at slurm file output, gives this error: terminate called after throwing an instance of 'std::out_of_range'
  what():  vector::_M_range_check: __n (which is 6259) >= this->size() (which is 6259)


```
nano star_align_5_2_test.sh
```

```
#!/bin/bash
#SBATCH --mem=64G

module load STAR/2.5.3a-foss-2016b


STAR --runMode alignReads --quantMode TranscriptomeSAM --outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/5_2_test_TMP --readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/5_2.fastq.trim.fq --genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_plob/ --sjdbGTFfeatureExon exon --sjdbGTFfile /data/putnamlab/jillashey/genome/Plutea/Plut.GFFannotation.fixed.gff --sjdbGTFtagExonParentTranscript ID --twopassMode Basic --twopass1readsN -1 --outStd Log BAM_Unsorted BAM_Quant --outSAMtype BAM Unsorted SortedByCoordinate --outReadsUnmapped Fastx --outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_plob/5_2_test.
```
Submitted batch job 1668262


Running without the gff and all other gff commands 

```
nano star_align_6_1_test_noGFF.sh

#!/bin/bash
#SBATCH --mem=64G

module load STAR/2.5.3a-foss-2016b


STAR --runMode alignReads --quantMode TranscriptomeSAM --outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/6_1_test_noGFF_TMP --readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/6_1.fastq.trim.fq --genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_plob/ --twopassMode Basic --twopass1readsN -1 --outStd Log BAM_Unsorted BAM_Quant --outSAMtype BAM Unsorted SortedByCoordinate --outReadsUnmapped Fastx --outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_plob/6_1_test_noGFF.
```
Submitted batch job 1668323

This worked! Going to try with another Plob sample to make sure it will really work

```
nano star_align_25_1_test_noGFF.sh

#!/bin/bash
#SBATCH --mem=64G

module load STAR/2.5.3a-foss-2016b


STAR --runMode alignReads --quantMode TranscriptomeSAM --outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/25_1_test_noGFF_TMP --readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/25_1.fastq.trim.fq --genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_plob/ --twopassMode Basic --twopass1readsN -1 --outStd Log BAM_Unsorted BAM_Quant --outSAMtype BAM Unsorted SortedByCoordinate --outReadsUnmapped Fastx --outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_plob/25_1_test_noGFF.
```
Submitted batch job 1668782

Yay worked for this one too. Now I'll submit an array job with all of the samples 


```
# Make directory for star alignment without GFF
mkdir staralign_plob_noGFF
cd staralign_plob_noGFF
nano star_align_plob_noGFF.sh
```

```
#!/bin/bash
#SBATCH --mem=64G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu

module load STAR/2.5.3a-foss-2016b

FILENAME=$(ls -l /data/putnamlab/jillashey/Francois_data/data/trimmed/ | sed -n $((${SLURM_ARRAY_TASK_ID}+1))p | awk '{print $9}' )

STAR --runMode alignReads --quantMode TranscriptomeSAM --outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/${FILENAME}_noGFF_TMP --readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/${FILENAME} --genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_plob/ --twopassMode Basic --twopass1readsN -1 --outStd Log BAM_Unsorted BAM_Quant --outSAMtype BAM Unsorted SortedByCoordinate --outReadsUnmapped Fastx --outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_plob/${FILENAME}.noGFF.
```

```
# Submit as an array job 
sbatch --array 1-64 star_align_plob_noGFF.sh
```
Submitted batch job 1668821

Weird it only says 33 of them are running. I'll see what happens in a few hrs
So only 33 of them ran. maybe thats the limit on array jobs?

Running the rest now 

```
# Submit as an array job 
sbatch --array 34-64 star_align_plob_noGFF.sh
```
Submitted batch job 1668912

Running some more that were not run
```
sbatch --array 62-64 star_align_plob_noGFF.sh
```
Submitted batch job 1668998

Still need 2, 26, 27, 28 in the order 

```
sbatch --array 2 star_align_plob_noGFF.sh
sbatch --array 26-28 star_align_plob_noGFF.sh
```

Submitted batch job 1669032
Submitted batch job 1669033

Running last one...not sure why they all got jumbled around

```
nano star_align_38_2_plut_noGFF.sh

#!/bin/bash
#SBATCH --mem=64G

module load STAR/2.5.3a-foss-2016b


STAR --runMode alignReads --quantMode TranscriptomeSAM --outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/38_2_noGFF_TMP --readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/38_2.fastq.trim.fq --genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_plob/ --twopassMode Basic --twopass1readsN -1 --outStd Log BAM_Unsorted BAM_Quant --outSAMtype BAM Unsorted SortedByCoordinate --outReadsUnmapped Fastx --outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_plob/38_2_noGFF.

sbatch star_align_38_2_plut_noGFF.sh
```

Submitted batch job 1669051

38_2 never finished, running again 

### xxxx) Perform gene counts with stringTie

#### P. acuta

a) Make folders for stringtie results 

```
cd /data/putnamlab/jillashey/Francois_data/Hawaii/stringTie_star
mkdir Pacuta
cd Pacuta
mkdir BAM GTF GTF_merge
```

b) Move BAM files to stringTie folder 

```
cd /data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/AlignReads_pacuta/pacuta_only
mv *Aligned.sortedByCoord.out.bam ../../../../stringTie_star/pacuta/BAM/
```

c) Assemble and estimate reads 

```
cd /data/putnamlab/jillashey/Francois_data/Hawaii/stringTie_star/pacuta/BAM

nano stringTie_pacuta_assemble.sh

#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=2
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="assemble_pacuta_out_error"
#SBATCH --output="assemble_pacuta_out"

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/Hawaii/stringTie_star/pacuta/BAM

array1=($(ls $F/*bam))
for i in ${array1[@]}; do
	stringtie -G /data/putnamlab/jillashey/genome/Pacuta/Pacuta.gff.annotations.fixed_transcript.gff3 -e -o ${i}.gtf ${i}
	echo "${i}"
done

sbatch stringTie_pacuta_assemble.sh
```

Submitted batch job 1960732

For some reason, sample 47_2 didn't run, no gtf file produced...making separate folder for 47_2 and running it by itself 

jk so for some reason the 47_2 BAM file is empty? going back to STAR to rerun 

Align 47_2 individually 

```
cd /data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/AlignReads_pacuta/pacuta_only

nano 47_2_align.sh

#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=2
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="47_2_align_pacuta_out_error"
#SBATCH --output="47_2_align_pacuta_out"

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads --quantMode TranscriptomeSAM --outTmpDir /data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/test_TMP --readFilesIn 47_2.fastq.trim.fq --genomeDir /data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/GenomeIndex_pacuta/ --twopassMode Basic --twopass1readsN -1 --outStd Log BAM_Unsorted BAM_Quant --outSAMtype BAM Unsorted SortedByCoordinate --outReadsUnmapped Fastx --outFileNamePrefix 47_2.

sbatch 47_2_align.sh
```

Submitted batch job 1960734

Move 47_2 BAM files to stringTie folder 

```
cd /data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/AlignReads_pacuta/pacuta_only
mv *Aligned.sortedByCoord.out.bam ../../../../stringTie_star/pacuta/BAM/
```

now assemble reads for sample 47_2. made its own folder to run separately 

```
cd /data/putnamlab/jillashey/Francois_data/Hawaii/stringTie_star/pacuta/BAM/47_2

nano stringTie_pacuta_assemble_47_2.sh

#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=2
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="assemble_pacuta_47_2_out_error"
#SBATCH --output="assemble_pacuta_47_2_out"

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/Hawaii/stringTie_star/pacuta/BAM/47_2

array1=($(ls $F/*bam))
for i in ${array1[@]}; do
        stringtie -G /data/putnamlab/jillashey/genome/Pacuta/Pacuta.gff.annotations.fixed_transcript.gff3 -e -o ${i}.gtf ${i}
        echo "${i}"
done

sbatch stringTie_pacuta_assemble_47_2.sh
```

Submitted batch job 1960750

Move 47_2 gtf file into GTF directory 

```
cd /data/putnamlab/jillashey/Francois_data/Hawaii/stringTie_star/pacuta/GTF

mv /data/putnamlab/jillashey/Francois_data/Hawaii/stringTie_star/pacuta/BAM/47_2/47_2.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf .
```

d) Merge stringTie gtf results 

```
ls *gtf > pacuta_mergelist.txt
cat pacuta_mergelist.txt

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

stringtie --merge -p 8 -G /data/putnamlab/jillashey/genome/Pacuta/Pacuta.gff.annotations.fixed_transcript.gff3 -o stringtie_pacuta_merged.gtf pacuta_mergelist.txt
```

e) Assess assembly quality

```
module load gffcompare/0.11.5-foss-2018b

gffcompare -r /data/putnamlab/jillashey/genome/Pacuta/Pacuta.gff.annotations.fixed_transcript.gff3 -o Pacuta.merged stringtie_pacuta_merged.gtf

38913 reference transcripts loaded.
  518 duplicate reference transcripts discarded.
  38672 query transfrags loaded.
```

f) Re-estimate assembly 

```
nano stringTie_pacuta_reassemble.sh

#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="reassemble_pacuta_out_error"
#SBATCH --output="reassemble_pacuta_out"

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/Hawaii/stringTie_star/pacuta/BAM

array1=($(ls $F/*bam))
for i in ${array1[@]}
do
stringtie -e -G /data/putnamlab/jillashey/genome/Pacuta/Pacuta.gff.annotations.fixed_transcript.gff3 -o ${i}.merge.gtf ${i}
echo "${i}"
done

sbatch stringTie_pacuta_reassemble.sh
```

Submitted batch job 1960752

Move merged gtf files 

```
mv *merge.gtf ../GTF_merge
```

g) Create gene matrix 


```
module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/Hawaii/stringTie_star/pacuta/GTF_merge/

array2=($(ls *merge.gtf))

for i in ${array2[@]}
do
echo "${i} $F${i}" >> sample_list_pacuta.txt
done

python prepDE.py -g gene_count_pacuta_matrix.csv -i sample_list_pacuta.txt

wc -l gene_count_pacuta_matrix.csv 
35709 gene_count_pacuta_matrix.csv
```

h) Secure-copy gene counts onto local computer

```
scp jillashey@bluewaves.uri.edu:/data/putnamlab/jillashey/Francois_data/Hawaii/stringTie_star/pacuta/GTF_merge/gene_count_pacuta_matrix.csv /Users/jillashey/Desktop/Putnamlab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/pacuta
```

The gene matrix has mostly STRG.xxx as gene ids (as opposed to the actual gene ids from the gff file). Going to try to use the [mstrg_prep.py](https://gist.github.com/gpertea/4207fa9cb30fe7fec0eb52bd29b9a976) python script to add gene ids to MSTRG gene ids. I think I need to run this on the stringtie_pacuta_merged.gtf file ??. See comments in this [thread](https://github.com/gpertea/stringtie/issues/113). 

This comment is of interest to me. So was I supposed to be using the stringtie_pacuta_merged.gtf for the next run of stringtie? 

```
gpertea commented on Aug 27, 2018

As explained in #179 (comment) , the idea was to use the script on the output of stringtie --merge, which is a GTF file (let's call it merged.gtf) and from it generate a "processed" gtf output (at stdout) which transforms the 'gene_id' attribute to include one or more reference gene_ids (or gene_name, if the script), if any were found to match any of the transcripts assembled there. There is a "usage" comment as the 2nd line of that gist suggesting the command line:

mstrg_prep.py merged.gtf > merged_prep.gtf
The input file is not directly modified -- a new gtf file is produced at stdout (captured as merged_prep.gtf in the example above), and that's the one that should be used with -G option for the subsequent runs of stringtie -e for each sample.
```

Downloaded the ```mstrg_prep.py``` script and going to run w/ merged gtf file 

```
module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

python mstrg_prep.py stringtie_pacuta_merged.gtf > merged_prep_test.gtf

```

Okay not quite sure what to do now. First, I'm going to try to re-estimate assembly with the stringtie_pacuta_merged.gtf file. 

```
nano stringTie_pacuta_reassemble_merge.sh

#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="reassemble_merge_pacuta_out_error"
#SBATCH --output="reassemble_merge_pacuta_out"

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/Hawaii/stringTie_star/pacuta/BAM

array1=($(ls $F/*bam))
for i in ${array1[@]}
do
stringtie -e -G /data/putnamlab/jillashey/Francois_data/Hawaii/stringTie_star/pacuta/GTF/stringtie_pacuta_merged.gtf -o ${i}.merge.gtf ${i}
echo "${i}"
done

sbatch stringTie_pacuta_reassemble_merge.sh
```

Submitted batch job 1960766

Move merged gtf files 

```
cd /data/putnamlab/jillashey/Francois_data/Hawaii/stringTie_star/pacuta/
mkdir GTF_merge_20220116
mv *merge.gtf ../GTF_merge_20220116
```

Create gene matrix

```
module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/Hawaii/stringTie_star/pacuta/GTF_merge_20220116/

array2=($(ls *merge.gtf))

for i in ${array2[@]}
do
echo "${i} $F${i}" >> sample_list_pacuta.txt
done

python prepDE.py -g gene_count_pacuta_matrix_merge.csv -i sample_list_pacuta.txt

wc -l 
```

Then, I'm going to re-estimate assembly using the merged_prep_test.gtf file and see how the gene count matrix files are different 

```
nano stringTie_pacuta_reassemble_merge_prep.sh

#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="reassemble_merge_prep_pacuta_out_error"
#SBATCH --output="reassemble_merge_prep_pacuta_out"

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/Hawaii/stringTie_star/pacuta/BAM

array1=($(ls $F/*bam))
for i in ${array1[@]}
do
stringtie -e -G /data/putnamlab/jillashey/Francois_data/Hawaii/stringTie_star/pacuta/GTF/merged_prep_test.gtf -o ${i}.merge_prep.gtf ${i}
echo "${i}"
done

sbatch stringTie_pacuta_reassemble_merge_prep.sh
```

Submitted batch job 1960768

Move merged gtf files 

```
mv *merge_prep.gtf ../GTF_merge_20220116
```

Create gene matrix

```
module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/Hawaii/stringTie_star/pacuta/GTF_merge_20220116/

array2=($(ls *merge_prep.gtf))

for i in ${array2[@]}
do
echo "${i} $F${i}" >> sample_list_prep_pacuta.txt
done

python prepDE.py -g gene_count_pacuta_matrix_merge_prep.csv -i sample_list_prep_pacuta.txt

wc -l 
```

Because of the number of MSTRG/STRG gene ids that are generated, I am going to use ```gene_count_pacuta_matrix_merge.csv``` to analyze gene counts. In R, I will merge this csv file w/ ```stringtie_pacuta_merged.gtf``` to match MSTRG/STRG ids w/ gene ids from original gff.

#### P. lobata

TAR was run without GFF file included 

a) Move reads (aligned by Coordinate) from output to stringTie folder 

For some reason, 38_2 not running. That's okay for now, as its a Pdam sample.

```
mv *Aligned.sortedByCoord.out.bam ../../../stringTie_star/plob/
```

b) Assemble reads with genome annotation

```
nano stringTie_plob_assemble.sh

#!/bin/bash
#SBATCH --mem=64GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --error="assemble_plob_out_error"
#SBATCH --output="assemble_plob_out"

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

# StringTie to assemble transcripts for each sample with the annotation file

F=/data/putnamlab/jillashey/Francois_data/Hawaii/stringTie_star/plob/BAM

array1=($(ls $F/*bam))
for i in ${array1[@]}; do
	stringtie -G /data/putnamlab/jillashey/genome/Plutea/Plut.GFFannotation.fixed_transcript.gff -e -o ${i}.gtf ${i}
	echo "${i}"
done

sbatch stringTie_plob_assemble.sh
```

Submitted batch job 1768981

c) Merge stringTie gtf results 

At this point, I'm only going to move forward with Plob samples

```
# move only plob samples
mv 6_1.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.gtf 6_2.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.gtf 7_1.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.gtf 8_1.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.gtf 9_1.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.gtf 9_2.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.gtf 21_1.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.gtf 22_1.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.gtf 23_1.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.gtf 25_1.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.gtf 25_2.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.gtf 26_1.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.gtf 26_2.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.gtf 27_1.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.gtf 29_1.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.gtf 34_1.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.gtf ../GTF 
```

```

ls *gtf > plob_mergelist.txt
cat plob_mergelist.txt

module load StringTie/2.1.1-GCCcore-7.3.0

stringtie --merge -p 8 -G /data/putnamlab/jillashey/genome/Plutea/Plut.GFFannotation.fixed_transcript.gff -o stringtie_plob_merged.gtf plob_mergelist.txt
```

d) Assess assembly quality

```
module load gffcompare/0.11.5-foss-2018b

gffcompare -r /data/putnamlab/jillashey/genome/Plutea/Plut.GFFannotation.fixed_transcript.gff -o Plob.merged stringtie_plob_merged.gtf

  31126 reference transcripts loaded.
  31126 query transfrags loaded.
```

e) Re-estimate assembly 

```
nano stringTie_plob_re-assemble.sh

#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="re-assemble_plob_out_error"
#SBATCH --output="re-assemble_plob_out"

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/Hawaii/stringTie_star/plob/BAM # do i need a / after BAM?

array1=($(ls $F/*bam))
for i in ${array1[@]}
do
stringtie -e -G /data/putnamlab/jillashey/genome/Plutea/Plut.GFFannotation.fixed_transcript.gff -o ${i}.merge.gtf ${i}
echo "${i}"
done

sbatch stringTie_plob_re-assemble.sh
```

Submitted batch job 1768989

```
# move only plob samples
mv 6_1.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.merge.gtf 6_2.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.merge.gtf 7_1.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.merge.gtf 8_1.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.merge.gtf 9_1.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.merge.gtf 9_2.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.merge.gtf 21_1.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.merge.gtf 22_1.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.merge.gtf 23_1.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.merge.gtf 25_1.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.merge.gtf 25_2.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.merge.gtf 26_1.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.merge.gtf 26_2.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.merge.gtf 27_1.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.merge.gtf 29_1.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.merge.gtf 34_1.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.merge.gtf ../GTF_merge 
```

f) Create gene matrix

```
module load StringTie/2.1.1-GCCcore-7.3.0
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/Hawaii/stringTie_star/plob/GTF_merge/

array2=($(ls *merge.gtf))

for i in ${array2[@]}
do
echo "${i} $F${i}" >> sample_list_plob.txt
done

python prepDE.py -g gene_count_plob_only_matrix.csv -i sample_list_plob.txt
```

g) Secure-copy gene counts onto local computer

```
scp jillashey@bluewaves.uri.edu:/data/putnamlab/jillashey/Francois_data/Hawaii/stringTie_star/plob/GTF_merge/gene_count_plob_only_matrix.csv /Users/jillashey/Desktop/Putnamlab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/
```