## STAR pacuta

For sediment stress paper, we started out analyzing the Pocillopora coral as P.dam. After further discussion, the Pocillopora coral is most likely P.acuta. So I need to re-run STAR using the Pacuta genomic info 

Using files from http://cyanophora.rutgers.edu/Pocillopora_acuta/. Has not been published yet. The files cannot be downloaded with 'wget' command, so the files need to be first downloaded onto local computer and then uploaded to server. 

File quality and trimming has already been done for these samples, so will go right into STAR

a) After downloading onto server, unzip files 

```
ls
braker_v1.codingseq.fasta.gz  braker_v1.gff3.gz  Pocillopora_acuta_HIv1.assembly.purged.fasta.gz  Pocillopora_acuta_HIv1.genes.pep.faa.gz

gunzip *
```

For some reason, some of the files are labelled 'braker', not sure what that means. I'm going to leave the names unchanged for now.

b) Generate genome index

Forgot that STAR needs transcript_id identifier to run properly, so added that identifier for Pacuta in R. Renamed gff file as Pacuta.gff.annotations.fixed_transcript.gff3

Create genome index file in ```/data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR
``` and run genome index step of STAR

```
mkdir GenomeIndex_pacuta

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

c) Align reads to genome 

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

Submitted batch job 1960686. cancelled job because it was pending for >1 day

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

d) Perform gene counts with stringTie

Make folders for stringtie results 

```
cd /data/putnamlab/jillashey/Francois_data/Hawaii/stringTie_star
mkdir Pacuta
cd Pacuta
mkdir BAM GTF GTF_merge
```

Move BAM files to stringTie folder 

```
cd /data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/AlignReads_pacuta/pacuta_only
mv *Aligned.sortedByCoord.out.bam ../../../../stringTie_star/pacuta/BAM/
```

Assemble and estimate reads 

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

Move gtf file into GTF directory 

```
cd /data/putnamlab/jillashey/Francois_data/Hawaii/stringTie_star/pacuta/GTF

mv /data/putnamlab/jillashey/Francois_data/Hawaii/stringTie_star/pacuta/BAM/47_2/47_2.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf .
```

Okay NOW I can merge stringTie gtf results 

```
ls *gtf > pacuta_mergelist.txt
cat pacuta_mergelist.txt

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

stringtie --merge -p 8 -G /data/putnamlab/jillashey/genome/Pacuta/Pacuta.gff.annotations.fixed_transcript.gff3 -o stringtie_pacuta_merged.gtf pacuta_mergelist.txt
```

Assess assembly quality

```
module load gffcompare/0.11.5-foss-2018b

gffcompare -r /data/putnamlab/jillashey/genome/Pacuta/Pacuta.gff.annotations.fixed_transcript.gff3 -o Pacuta.merged stringtie_pacuta_merged.gtf

38913 reference transcripts loaded.
  518 duplicate reference transcripts discarded.
  38672 query transfrags loaded.
```

Re-estimate assembly 

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

Create gene matrix

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

Secure-copy gene counts onto local computer

```
scp jillashey@bluewaves.uri.edu:/data/putnamlab/jillashey/Francois_data/Hawaii/stringTie_star/pacuta/GTF_merge/gene_count_pacuta_matrix.csv /Users/jillashey/Desktop/Putnamlab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/pacuta
```