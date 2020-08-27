
a) make directories for STAR output
```
```

b) Obtain reference genome (.fasta or .fna) and reference annotation (.gtf or .gff). Unzip these files

```
# P.dam genome info
/data/putnamlab/REFS/Pdam
# GCF_003704095.1_ASM370409v1_genomic.fna.gz -- genome file
# GCF_003704095.1_ASM370409v1_genomic.gff.gz -- annotation file
```
For some reason, REFS won't let me access the genome info in that folder. I moved them into my own genome folder so I can access them

```
# Pdam unzip
gunzip GCF_003704095.1_ASM370409v1_genomic.fna.gz
gunzip GCF_003704095.1_ASM370409v1_genomic.gff
```

c) Unzip fastq files so they can properly analyzed by STAR. Sort files by species 

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

```
module load STAR/2.5.3a-foss-2016b

STAR --runThreadN 10 \
--runMode genomeGenerate \
-- genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pdam/ \
--genomeFastaFiles /data/putnamlab/jillashey/genome/Pdam/GCF_003704095.1_ASM370409v1_genomic.fna \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pdam/GCF_003704095.1_ASM370409v1_genomic.gff

Jul 29 22:31:38 ..... started STAR run
Jul 29 22:31:38 ... starting to generate Genome files
Jul 29 22:33:16 ... starting to sort Suffix Array. This may take a long time...
Jul 29 22:33:19 ... sorting Suffix Array chunks and saving them to disk...
Jul 29 22:34:48 ... loading chunks from disk, packing SA...
Jul 29 22:36:11 ... finished generating suffix array
Jul 29 22:36:11 ... generating Suffix Array index
Jul 29 22:37:48 ... completed Suffix Array index
Jul 29 22:37:48 ..... processing annotations GTF
Jul 29 22:38:16 ..... inserting junctions into the genome indices
Jul 29 22:39:56 ... writing Genome to disk ...
Jul 29 22:40:04 ... writing Suffix Array to disk ...
Jul 29 22:40:20 ... writing SAindex to disk
Jul 29 22:40:31 ..... finished successfully
```

e) Align reads 

For P.dam

Attempt 1

```
nano star_align_batch_script_pdam.sh
```

```
#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

FILENAME=$(ls -l /data/putnamlab/jillashey/Francois_data/data/trimmed/Pdam/ | sed -n $((${SLURM_ARRAY_TASK_ID}+1))p | awk '{print $9}' )

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/${FILENAME}_TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/Pdam/${FILENAME} \ 
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pdam/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript Parent \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pdam/GCF_003704095.1_ASM370409v1_genomic.gff
--twopassMode Basic \ 
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \ 
--outSAMtype BAM Unsorted \ # add sorted here?
--outReadsUnmapped Fastx \ 
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pdam/${FILENAME}

```

```
# Submit as an array job 
sbatch --array star_align_batch_script_pdam.sh
```

```
EXITING because of fatal ERROR: could not make temporary directory: /data/putnamlab/jillashey/Francois_data/output/STAR/
SOLUTION: (i) please check the path and writing permissions 
 (ii) if you specified --outTmpDir, and this directory exists - please remove it before running STAR

Jul 31 14:16:52 ...... FATAL ERROR, exiting
/var/spool/slurmd/job1663148/slurm_script: line 11: --genomeDir: command not found
/var/spool/slurmd/job1663148/slurm_script: line 15: --twopassMode: command not found
/var/spool/slurmd/job1663148/slurm_script: line 16: --twopass1readsN: command not found
/var/spool/slurmd/job1663148/slurm_script: line 18: --outSAMtype: command not found
/var/spool/slurmd/job1663148/slurm_script: line 19: --outReadsUnmapped: command not found
/var/spool/slurmd/job1663148/slurm_script: line 20: --outFileNamePrefix: command not found
```

Going to try without --tmpDir

Attempt 2

```
nano star_align_no_tmpDir_pdam.sh

#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

FILENAME=$(ls -l /data/putnamlab/jillashey/Francois_data/data/trimmed/Pdam/ | sed -n $((${SLURM_ARRAY_TASK_ID}+1))p | awk '{print $9}' )

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/Pdam/${FILENAME} \ 
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pdam/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript Parent \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pdam/GCF_003704095.1_ASM370409v1_genomic.gff
--twopassMode Basic \ 
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \ 
--outSAMtype BAM Unsorted \ # add sorted here?
--outReadsUnmapped Fastx \ 
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pdam/${FILENAME}

sbatch --array 1-2 star_align_batch_script_pdam.sh
```

When I do that, gives me error: cannot read files in

Going to try with only one sample now 

Attempt 3

```
nano star_align_one_sample.sh

#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads --quantMode TranscriptomeSAM --readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/Pdam/10_2.fastq.trim.fq --genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pdam/ --sjdbGTFfeatureExon exon --sjdbGTFtagExonParentTranscript Parent --sjdbGTFfile /data/putnamlab/jillashey/genome/Pdam/GCF_003704095.1_ASM370409v1_genomic.gff --twopassMode Basic --twopass1readsN -1 --outStd Log BAM_Unsorted BAM_Quant --outSAMtype BAM Unsorted --outReadsUnmapped Fastx --outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pdam/10_2.fastq.trim.fq

sbatch star_align_one_sample.sh
```

Stuck here: 
Jul 31 15:16:08 ..... started STAR run
Jul 31 15:16:08 ..... loading genome
Jul 31 15:16:11 ..... processing annotations GTF
Jul 31 15:16:38 ..... inserting junctions into the genome indices
Jul 31 15:16:51 ..... started 1st pass mapping
client_loop: send disconnect: Broken pipe

Will not move past 1st pass mapping. Not sure why because the pdam genome/annotations worked fine with connelly data.

Potential solutions (from alex love github): 

- the problem is not in smallness of the genome per se, but rather in the incompleteness of the genome. Since majority of the reads will not map to the rRNA reference, STAR would be trying hard to place them, which slows down the mapping speed.
I think including the rRNA sequences with the genome is the best option. It will increase the multi-mapping for some reads, however, can deal with such alignments in the postprocessing. Most of rRNA alignments will be multi-mappers anyway.
Another option is to map to the standard genome first, and then re-map the unmapped reads (output with --outReadsUnmapped Fastx) to the rRNA reference. --> though I dont't think this is it because the pdam genome is pretty good
- Large number of junctions (>500k) will slow down the 2nd pass mapping. Another possibility that there is a highly expressed locus with many novel junctions, typically false positive junctions - e.g. in mitochondrion genome, or in rRNA loci.
	- In any case, you would need to filter the SJ.out.tab file before the 2nd pass mapping, so you would need to do "manual" 2-pass mapping (without --twopassMode option):

1)Map 1st pass with standard parameters, without BAM/SAM output (--outSAMtype None), without --chim* options and without --outFilterType BySJout.
2)Filter the resulting SJ.out.tab file. I would start with removing all junctions mapping to chrM and non-chromosomal contigs.
3)Run second pass (in a separate directory) with all the output options and --sjdbFileChrStartEnd /path/to/SJ.out.tab.filtered.

Attempt 4

```
nano star_align_batch_script_pdam.sh
```

```
#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

FILENAME=$(ls -l /data/putnamlab/jillashey/Francois_data/data/trimmed/ | sed -n $((${SLURM_ARRAY_TASK_ID}+1))p | awk '{print $9}' )

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/${FILENAME}_TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/${FILENAME} \ 
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pdam/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript Parent \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pdam/GCF_003704095.1_ASM370409v1_genomic.gff \
--sjdbGTFtagExonParentTranscript ID \
--twopassMode Basic \ 
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \ 
--outSAMtype BAM Unsorted SortedByCoordinate\ 
--outReadsUnmapped Fastx \ 
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pdam/${FILENAME}.

```

```
# Submit as an array job 
sbatch --array 1-32 star_align_batch_script_pdam.sh
```
Submitted batch job 1663852

Attempt 5

Need to do it with all 64 samples

```
nano star_align_batch_script_pdam.sh
```

```
#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

FILENAME=$(ls -l /data/putnamlab/jillashey/Francois_data/data/trimmed/ | sed -n $((${SLURM_ARRAY_TASK_ID}+1))p | awk '{print $9}' )

echo $FILENAME

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/${FILENAME}_TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/${FILENAME} \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pdam/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pdam/GCF_003704095.1_ASM370409v1_genomic.gff \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pdam/${FILENAME}.
```

```
sbatch --array 1-64 star_align_batch_script_pdam.sh
```

Submitted batch job 1667070

Worked for the majority of samples, but some were not run for some reason. Just going to run them separately because its not too many

```
nano star_23_1_pdam.sh

#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/23_1.fastq.trim.fq _TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/23_1.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pdam/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pdam/GCF_003704095.1_ASM370409v1_genomic.gff \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pdam/23_1.fastq.trim.fq 

sbatch star_23_1_pdam.sh
```
Submitted batch job 1667276


```
nano star_24_1_pdam.sh

#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/24_1.fastq.trim.fq _TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/24_1.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pdam/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pdam/GCF_003704095.1_ASM370409v1_genomic.gff \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pdam/24_1.fastq.trim.fq 

sbatch star_24_1_pdam.sh
```

Submitted batch job 1667287


```
nano star_25_1_pdam.sh

#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/25_1.fastq.trim.fq _TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/25_1.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pdam/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pdam/GCF_003704095.1_ASM370409v1_genomic.gff \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pdam/25_1.fastq.trim.fq 

sbatch star_25_1_pdam.sh
```

Submitted batch job 1667288



```
nano star_26_1_pdam.sh

#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/26_1.fastq.trim.fq _TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/26_1.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pdam/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pdam/GCF_003704095.1_ASM370409v1_genomic.gff \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pdam/26_1.fastq.trim.fq 

sbatch star_26_1_pdam.sh
```

Submitted batch job 1667289



```
nano star_27_1_pdam.sh

#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/27_1.fastq.trim.fq _TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/27_1.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pdam/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pdam/GCF_003704095.1_ASM370409v1_genomic.gff \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pdam/27_1.fastq.trim.fq 

sbatch star_27_1_pdam.sh
```

Submitted batch job 1667290



```
nano star_28_2_pdam.sh

#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/28_2.fastq.trim.fq _TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/28_2.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pdam/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pdam/GCF_003704095.1_ASM370409v1_genomic.gff \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pdam/28_2.fastq.trim.fq 

sbatch star_28_2_pdam.sh
```

Submitted batch job 1667291



```
nano star_29_1_pdam.sh

#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/29_1.fastq.trim.fq _TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/29_1.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pdam/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pdam/GCF_003704095.1_ASM370409v1_genomic.gff \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pdam/29_1.fastq.trim.fq 

sbatch star_29_1_pdam.sh
```

Submitted batch job 1667292



```
nano star_30_1_pdam.sh

#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/30_1.fastq.trim.fq _TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/30_1.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pdam/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pdam/GCF_003704095.1_ASM370409v1_genomic.gff \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pdam/30_1.fastq.trim.fq 

sbatch star_30_1_pdam.sh
```

Submitted batch job 1667293



```
nano star_31_2_pdam.sh

#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/31_2.fastq.trim.fq _TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/31_2.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pdam/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pdam/GCF_003704095.1_ASM370409v1_genomic.gff \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pdam/31_2.fastq.trim.fq 

sbatch star_31_2_pdam.sh
```

Submitted batch job 1667294

```
nano star_25_2_pdam.sh

#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/25_2.fastq.trim.fq _TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/25_2.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pdam/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pdam/GCF_003704095.1_ASM370409v1_genomic.gff \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pdam/25_2.fastq.trim.fq 

sbatch star_25_2_pdam.sh
```
Submitted batch job 1667444


```
nano star_26_2_pdam.sh

#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/26_2.fastq.trim.fq _TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/26_2.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pdam/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pdam/GCF_003704095.1_ASM370409v1_genomic.gff \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pdam/26_2.fastq.trim.fq 

sbatch star_26_2_pdam.sh
```
Submitted batch job 1667445

Going to try to run a pdam sample without the GFF files to see how it compares 

```
nano star_align_pdam_1_2_test_noGFF.sh

#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/1_2_test_noGFF_TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/1_2.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pdam/ \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pdam/1_2_test_noGFF.

sbatch star_align_pdam_1_2_test_noGFF.sh
```

Submitted batch job 1669003

Same output as STAR alignment with the GTF file, good sign
