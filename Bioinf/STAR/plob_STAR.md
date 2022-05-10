# P. lobata (using P. lutea genome)

a) make directories for STAR output

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

STAR --runThreadN 10 --runMode genomeGenerate --genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_plob/ --genomeFastaFiles /data/putnamlab/jillashey/genome/Plutea/plut_final_2.1.fasta --sjdbGTFfile /data/putnamlab/jillashey/genome/Plutea/Plut.GFFannotation.fixed_transcript.gff
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

StringTie - STAR, plobata samples with plutea genome 

STAR was run without GFF file included 

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

20220116

Going to try to run re-assemble step with the stringtie_plob_merged.gtf file instead of the Plut.GFFannotation.fixed_transcript.gff file. 

```
cd /data/putnamlab/jillashey/Francois_data/Hawaii/stringTie_star/plob/BAM/plob_only

nano stringTie_plob_re-assemble_merge.sh

#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="re-assemble_merge_plob_out_error"
#SBATCH --output="re-assemble_merge_plob_out"

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/Hawaii/stringTie_star/plob/BAM/plob_only # do i need a / after BAM?

array1=($(ls $F/*bam))
for i in ${array1[@]}
do
stringtie -e -G /data/putnamlab/jillashey/Francois_data/Hawaii/stringTie_star/plob/GTF/stringtie_plob_merged.gtf -o ${i}.merge.gtf ${i}
echo "${i}"
done

sbatch stringTie_plob_re-assemble_merge.sh
```

Submitted batch job 1960767

Move samples to new GTF folder 

f) Create gene matrix

```
module load StringTie/2.1.1-GCCcore-7.3.0
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/Hawaii/stringTie_star/plob/GTF_merge_20220116/

array2=($(ls *merge.gtf))

for i in ${array2[@]}
do
echo "${i} $F${i}" >> sample_list_plob.txt
done

python prepDE.py -g gene_count_plob_matrix_merge.csv -i sample_list_plob.txt

```

I'm sticking w/ the original gene count matrix that I made 
