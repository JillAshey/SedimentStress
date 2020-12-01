# P. compressa

a-c) in other files but will add here later

d) Generate genome index

In order to have the gff file work with STAR, had to convert the ID in 'gene' column to transcript_id. R script for Pcomp ID conversion is here. 

Attempt 1

```
# getting new gff file from local computer
scp /Users/jillashey/Desktop/Pcomp.GFFannotation.fixed.gff jillashey@bluewaves.uri.edu:/data/putnamlab/jillashey/genome/Pcomp/
```

```
module load STAR/2.5.3a-foss-2016b

STAR --runThreadN 10 --runMode genomeGenerate --genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pcomp/ --genomeFastaFiles /data/putnamlab/jillashey/genome/Pcomp/Porites_compressa_contigs.fasta --sjdbGTFfile /data/putnamlab/jillashey/genome/Pcomp/Pcomp.GFFannotation.fixed.gff
```

Got this error:
Fatal INPUT FILE error, no exon lines in the GTF file: /data/putnamlab/jillashey/genome/Mcap/Mcap.GFFannotation.gff
Solution: check the formatting of the GTF file, it must contain some lines with exon in the 3rd column.
          Make sure the GTF file is unzipped.
          If exons are marked with a different word, use --sjdbGTFfeatureExon .

Error happening because the features (or id) column of Mcap annotation only has CSD, introns, and genes in it. STAR needs exon features to run. 
Possible solution: if you have no splices (i.e. all transcript are single-exon), you do not need to use the annotations (GTF file) at all. Or, indeed, you can use a non-splice-aware aligner like bwa or bowtie.
If you want to count reads per gene, you can rename all features in your GTF as "exon".

Attempt 2
STAR needs transcript ID to generate genome index. In R, we removed some info from ID and added it to transcript id, which was added to gene column

```
# getting new gff file from local computer
scp /Users/jillashey/Desktop/Pcomp.GFFannotation.fixed_transcript.gff jillashey@bluewaves.uri.edu:/data/putnamlab/jillashey/genome/Pcomp/
```

```
module load STAR/2.5.3a-foss-2016b

STAR --runThreadN 10 --runMode genomeGenerate --genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pcomp/ --genomeFastaFiles /data/putnamlab/jillashey/genome/Pcomp/Porites_compressa_contigs.fasta --sjdbGTFfile /data/putnamlab/jillashey/genome/Pcomp/Pcomp.GFFannotation.fixed_transcript.gff
```

:%s#old#new#g -- for vim replacing words 

e) Align reads 

Try with one sample first 

```
nano star_align_pcomp_test.sh
```

```
#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/test_14_2_TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/14_2.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pcomp/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pcomp/Pcomp.GFFannotation.fixed_transcript.gff
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pcomp/test_14_2.
```

```
sbatch star_align_pcomp_test.sh
```

Submitted batch job 1667423

Worked! Got a log.final.out.file, so will move forward will all 64 samples 

Attempt 1)

```
nano star_align_batch_script_pcomp.sh
```

```
#!/bin/bash
#SBATCH --mem=64GB

module load STAR/2.5.3a-foss-2016b

FILENAME=$(ls -l /data/putnamlab/jillashey/Francois_data/data/trimmed/ | sed -n $((${SLURM_ARRAY_TASK_ID}+1))p | awk '{print $9}' )

echo $FILENAME

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/${FILENAME}_TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/${FILENAME} \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pcomp/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pcomp/Pcomp.GFFannotation.fixed_transcript.gff
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pcomp/${FILENAME}.
```

```
sbatch --array 1-64 star_align_batch_script_pcomp.sh
```

Submitted batch job 1667451

Most were not processed because of error: EXITING because of fatal ERROR: could not make run-time genome directory directory: ./_STARgenome/  

Attempt 2)

Going to put all arguments together instead of on their own lines. Sometimes its just a spacing issue

```
#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

FILENAME=$(ls -l /data/putnamlab/jillashey/Francois_data/data/trimmed/ | sed -n $((${SLURM_ARRAY_TASK_ID}+1))p | awk '{print $9}' )

echo $FILENAME

STAR --runMode alignReads --quantMode TranscriptomeSAM --outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/${FILENAME}_TMP --readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/${FILENAME} --genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pcomp/ --sjdbGTFfeatureExon exon --sjdbGTFtagExonParentTranscript ID --sjdbGTFfile /data/putnamlab/jillashey/genome/Pcomp/Pcomp.GFFannotation.fixed_transcript.gff --twopassMode Basic --twopass1readsN -1 --outStd Log BAM_Unsorted BAM_Quant --outSAMtype BAM Unsorted SortedByCoordinate --outReadsUnmapped Fastx --outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pcomp/${FILENAME}.
```

```
sbatch --array 1-64 star_align_batch_script_pcomp.sh
```

Submitted batch job 1667546

Ones that need to be run individually

1) 

```
nano star_1_2_pcomp.sh

#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/1_2.fastq.trim.fq _TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/1_2.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pcomp/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pcomp/Pcomp.GFFannotation.fixed_transcript.gff \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pcomp/IndividualStarRuns/1_2.fastq.trim.fq 

sbatch star_1_2_pcomp.sh
```
Submitted batch job 1667641

2)

```
nano star_2_2_pcomp.sh

#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/2_2.fastq.trim.fq _TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/2_2.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pcomp/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pcomp/Pcomp.GFFannotation.fixed_transcript.gff \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pcomp/IndividualStarRuns/2_2.fastq.trim.fq 

sbatch star_2_2_pcomp.sh
```

Submitted batch job 1667642

3)

```
nano star_3_1_pcomp.sh

#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/3_1.fastq.trim.fq _TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/3_1.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pcomp/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pcomp/Pcomp.GFFannotation.fixed_transcript.gff \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pcomp/IndividualStarRuns/3_1.fastq.trim.fq 

sbatch star_3_1_pcomp.sh
```
Submitted batch job 1667643


4)

```
nano star_3_2_pcomp.sh

#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/3_2.fastq.trim.fq _TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/3_2.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pcomp/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pcomp/Pcomp.GFFannotation.fixed_transcript.gff \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pcomp/IndividualStarRuns/3_2.fastq.trim.fq 

sbatch star_3_2_pcomp.sh
```
Submitted batch job 1667644


5)

```
nano star_4_1_pcomp.sh

#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/4_1.fastq.trim.fq _TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/4_1.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pcomp/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pcomp/Pcomp.GFFannotation.fixed_transcript.gff \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pcomp/IndividualStarRuns/4_1.fastq.trim.fq 

sbatch star_4_1_pcomp.sh
```

Submitted batch job 1667645

6)

```
nano star_5_1_pcomp.sh

#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/5_1.fastq.trim.fq _TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/5_1.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pcomp/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pcomp/Pcomp.GFFannotation.fixed_transcript.gff \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pcomp/IndividualStarRuns/5_1.fastq.trim.fq 

sbatch star_5_1_pcomp.sh
```
Submitted batch job 1667646


7)

```
nano star_5_2_pcomp.sh

#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/5_2.fastq.trim.fq _TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/5_2.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pcomp/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pcomp/Pcomp.GFFannotation.fixed_transcript.gff \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pcomp/IndividualStarRuns/5_2.fastq.trim.fq 

sbatch star_5_2_pcomp.sh
```
Submitted batch job 1667647


8)

```
nano star_6_1_pcomp.sh

#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/6_1.fastq.trim.fq _TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/6_1.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pcomp/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pcomp/Pcomp.GFFannotation.fixed_transcript.gff \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pcomp/IndividualStarRuns/6_1.fastq.trim.fq 

sbatch star_6_1_pcomp.sh
```
Submitted batch job 1667648


9) 

```
nano star_7_1_pcomp.sh

#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/7_1.fastq.trim.fq _TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/7_1.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pcomp/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pcomp/Pcomp.GFFannotation.fixed_transcript.gff \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pcomp/IndividualStarRuns/7_1.fastq.trim.fq 

sbatch star_7_1_pcomp.sh
```

Submitted batch job 1667649


10)

```
nano star_9_1_pcomp.sh

#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/9_1.fastq.trim.fq _TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/9_1.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pcomp/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pcomp/Pcomp.GFFannotation.fixed_transcript.gff \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pcomp/IndividualStarRuns/9_1.fastq.trim.fq 

sbatch star_9_1_pcomp.sh
```

Submitted batch job 1667650


11)

```
nano star_9_2_pcomp.sh

#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/9_2.fastq.trim.fq _TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/9_2.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pcomp/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pcomp/Pcomp.GFFannotation.fixed_transcript.gff \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pcomp/IndividualStarRuns/9_2.fastq.trim.fq 

sbatch star_9_2_pcomp.sh
```
Submitted batch job 1667651


12) 

```
nano star_10_2_pcomp.sh

#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/10_2.fastq.trim.fq _TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/10_2.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pcomp/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pcomp/Pcomp.GFFannotation.fixed_transcript.gff \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pcomp/IndividualStarRuns/10_2.fastq.trim.fq 

sbatch star_10_2_pcomp.sh
```
Submitted batch job 1667652

13) 

```
nano star_11_1_pcomp.sh

#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/11_1.fastq.trim.fq _TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/11_1.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pcomp/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pcomp/Pcomp.GFFannotation.fixed_transcript.gff \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pcomp/IndividualStarRuns/11_1.fastq.trim.fq 

sbatch star_11_1_pcomp.sh
```
Submitted batch job 1667656


14) 

```
nano star_11_2_pcomp.sh

#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/11_2.fastq.trim.fq _TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/11_2.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pcomp/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pcomp/Pcomp.GFFannotation.fixed_transcript.gff \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pcomp/IndividualStarRuns/11_2.fastq.trim.fq 

sbatch star_11_2_pcomp.sh
```

Submitted batch job 1667657


15)

```
nano star_13_1_pcomp.sh

#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/13_1.fastq.trim.fq _TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/13_1.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pcomp/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pcomp/Pcomp.GFFannotation.fixed_transcript.gff \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pcomp/IndividualStarRuns/13_1.fastq.trim.fq 

sbatch star_13_1_pcomp.sh
```
Submitted batch job 1667658


16)

```
nano star_17_1_pcomp.sh

#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/17_1.fastq.trim.fq _TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/17_1.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pcomp/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pcomp/Pcomp.GFFannotation.fixed_transcript.gff \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pcomp/IndividualStarRuns/17_1.fastq.trim.fq 

sbatch star_17_1_pcomp.sh
```

Submitted batch job 1667659


17)


```
nano star_19_2_pcomp.sh

#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/19_2.fastq.trim.fq _TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/19_2.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pcomp/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pcomp/Pcomp.GFFannotation.fixed_transcript.gff \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pcomp/IndividualStarRuns/19_2.fastq.trim.fq 

sbatch star_19_2_pcomp.sh
```
Submitted batch job 1667660


18)

```
nano star_20_2_pcomp.sh

#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/20_2.fastq.trim.fq _TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/20_2.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pcomp/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pcomp/Pcomp.GFFannotation.fixed_transcript.gff \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pcomp/IndividualStarRuns/20_2.fastq.trim.fq 

sbatch star_20_2_pcomp.sh
```

Submitted batch job 1667661


19)

```
nano star_21_1_pcomp.sh

#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/21_1.fastq.trim.fq _TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/21_1.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pcomp/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pcomp/Pcomp.GFFannotation.fixed_transcript.gff \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pcomp/IndividualStarRuns/21_1.fastq.trim.fq 

sbatch star_21_1_pcomp.sh
```

Submitted batch job 1667662


20) 

```
nano star_22_1_pcomp.sh

#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/22_1.fastq.trim.fq _TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/22_1.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pcomp/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pcomp/Pcomp.GFFannotation.fixed_transcript.gff \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pcomp/IndividualStarRuns/22_1.fastq.trim.fq 

sbatch star_22_1_pcomp.sh
```

Submitted batch job 1667663


21) 

```
nano star_23_1_pcomp.sh

#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/23_1.fastq.trim.fq _TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/23_1.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pcomp/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pcomp/Pcomp.GFFannotation.fixed_transcript.gff \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pcomp/IndividualStarRuns/23_1.fastq.trim.fq 

sbatch star_23_1_pcomp.sh
```

Submitted batch job 1667664


22) 

```
nano star_24_1_pcomp.sh

#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/24_1.fastq.trim.fq _TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/24_1.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pcomp/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pcomp/Pcomp.GFFannotation.fixed_transcript.gff \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pcomp/IndividualStarRuns/24_1.fastq.trim.fq 

sbatch star_24_1_pcomp.sh
```

Submitted batch job 1667665


23) 

```
nano star_25_1_pcomp.sh

#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/25_1.fastq.trim.fq _TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/25_1.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pcomp/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pcomp/Pcomp.GFFannotation.fixed_transcript.gff \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pcomp/IndividualStarRuns/25_1.fastq.trim.fq 

sbatch star_25_1_pcomp.sh
```

Submitted batch job 1667666


24)

```
nano star_26_1_pcomp.sh

#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/26_1.fastq.trim.fq _TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/26_1.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pcomp/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pcomp/Pcomp.GFFannotation.fixed_transcript.gff \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pcomp/IndividualStarRuns/26_1.fastq.trim.fq 

sbatch star_26_1_pcomp.sh
```
Submitted batch job 1667667


25)

```
nano star_27_1_pcomp.sh

#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/27_1.fastq.trim.fq _TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/27_1.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pcomp/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pcomp/Pcomp.GFFannotation.fixed_transcript.gff \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pcomp/IndividualStarRuns/27_1.fastq.trim.fq 

sbatch star_27_1_pcomp.sh
```

Submitted batch job 1667668


26)

```
nano star_28_2_pcomp.sh

#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/28_2.fastq.trim.fq _TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/28_2.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pcomp/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pcomp/Pcomp.GFFannotation.fixed_transcript.gff \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pcomp/IndividualStarRuns/28_2.fastq.trim.fq 

sbatch star_28_2_pcomp.sh
```

Submitted batch job 1667669

27) 

```
nano star_29_1_pcomp.sh

#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/29_1.fastq.trim.fq _TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/29_1.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pcomp/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pcomp/Pcomp.GFFannotation.fixed_transcript.gff \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pcomp/IndividualStarRuns/29_1.fastq.trim.fq 

sbatch star_29_1_pcomp.sh
```

Submitted batch job 1667672

28) 

```
nano star_30_1_pcomp.sh

#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/30_1.fastq.trim.fq _TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/30_1.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pcomp/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pcomp/Pcomp.GFFannotation.fixed_transcript.gff \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pcomp/IndividualStarRuns/30_1.fastq.trim.fq 

sbatch star_30_1_pcomp.sh
```

Submitted batch job 1667673

29)

```
nano star_31_2_pcomp.sh

#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/31_2.fastq.trim.fq _TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/31_2.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pcomp/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pcomp/Pcomp.GFFannotation.fixed_transcript.gff \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pcomp/IndividualStarRuns/31_2.fastq.trim.fq 

sbatch star_31_2_pcomp.sh
```

Submitted batch job 1667674

30)

```
nano star_32_1_pcomp.sh

#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/32_1.fastq.trim.fq _TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/32_1.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pcomp/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pcomp/Pcomp.GFFannotation.fixed_transcript.gff \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pcomp/IndividualStarRuns/32_1.fastq.trim.fq 

sbatch star_32_1_pcomp.sh
```

Submitted batch job 1667675

31)

```
nano star_32_2_pcomp.sh

#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/32_2.fastq.trim.fq _TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/32_2.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pcomp/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pcomp/Pcomp.GFFannotation.fixed_transcript.gff \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pcomp/IndividualStarRuns/32_2.fastq.trim.fq 

sbatch star_32_2_pcomp.sh
```

Submitted batch job 1667676

32)

```
nano star_34_1_pcomp.sh

#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/34_1.fastq.trim.fq _TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/34_1.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pcomp/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pcomp/Pcomp.GFFannotation.fixed_transcript.gff \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pcomp/IndividualStarRuns/34_1.fastq.trim.fq 

sbatch star_34_1_pcomp.sh
```
Submitted batch job 1667677


33)

```
nano star_35_2_pcomp.sh

#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/35_2.fastq.trim.fq _TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/35_2.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pcomp/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pcomp/Pcomp.GFFannotation.fixed_transcript.gff \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pcomp/IndividualStarRuns/35_2.fastq.trim.fq 

sbatch star_35_2_pcomp.sh
```

Submitted batch job 1667679


34)


```
nano star_37_1_pcomp.sh

#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/37_1.fastq.trim.fq _TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/37_1.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pcomp/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pcomp/Pcomp.GFFannotation.fixed_transcript.gff \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pcomp/IndividualStarRuns/37_1.fastq.trim.fq 

sbatch star_37_1_pcomp.sh
```
Submitted batch job 1667680


35)

```
nano star_37_2_pcomp.sh

#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/37_2.fastq.trim.fq _TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/37_2.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pcomp/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pcomp/Pcomp.GFFannotation.fixed_transcript.gff \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pcomp/IndividualStarRuns/37_2.fastq.trim.fq 

sbatch star_37_2_pcomp.sh
```
Submitted batch job 1667681


36)


```
nano star_39_1_pcomp.sh

#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/39_1.fastq.trim.fq _TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/39_1.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pcomp/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pcomp/Pcomp.GFFannotation.fixed_transcript.gff \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pcomp/IndividualStarRuns/39_1.fastq.trim.fq 

sbatch star_39_1_pcomp.sh
```
Submitted batch job 1667682


37)

```
nano star_39_2_pcomp.sh

#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/39_2.fastq.trim.fq _TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/39_2.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pcomp/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pcomp/Pcomp.GFFannotation.fixed_transcript.gff \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pcomp/IndividualStarRuns/39_2.fastq.trim.fq 

sbatch star_39_2_pcomp.sh
```
Submitted batch job 1667683


38)

```
nano star_41_1_pcomp.sh

#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/41_1.fastq.trim.fq _TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/41_1.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pcomp/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pcomp/Pcomp.GFFannotation.fixed_transcript.gff \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pcomp/IndividualStarRuns/41_1.fastq.trim.fq 

sbatch star_41_1_pcomp.sh
```
Submitted batch job 1667684


39)

```
nano star_41_2_pcomp.sh

#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/41_2.fastq.trim.fq _TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/41_2.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pcomp/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pcomp/Pcomp.GFFannotation.fixed_transcript.gff \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pcomp/IndividualStarRuns/41_2.fastq.trim.fq 

sbatch star_41_2_pcomp.sh
```
Submitted batch job 1667685


40)

```
nano star_42_2_pcomp.sh

#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/42_2.fastq.trim.fq _TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/42_2.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pcomp/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pcomp/Pcomp.GFFannotation.fixed_transcript.gff \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pcomp/IndividualStarRuns/42_2.fastq.trim.fq 

sbatch star_42_2_pcomp.sh
```

Submitted batch job 1667686


41)

```
nano star_43_1_pcomp.sh

#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/43_1.fastq.trim.fq _TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/43_1.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pcomp/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pcomp/Pcomp.GFFannotation.fixed_transcript.gff \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pcomp/IndividualStarRuns/43_1.fastq.trim.fq 

sbatch star_43_1_pcomp.sh
```
Submitted batch job 1667687


42)

```
nano star_43_2_pcomp.sh

#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/43_2.fastq.trim.fq _TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/43_2.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pcomp/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pcomp/Pcomp.GFFannotation.fixed_transcript.gff \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pcomp/IndividualStarRuns/43_2.fastq.trim.fq 

sbatch star_43_2_pcomp.sh
```
Submitted batch job 1667689


43)

```
nano star_45_1_pcomp.sh

#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/45_1.fastq.trim.fq _TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/45_1.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pcomp/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pcomp/Pcomp.GFFannotation.fixed_transcript.gff \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pcomp/IndividualStarRuns/45_1.fastq.trim.fq 

sbatch star_45_1_pcomp.sh
```
Submitted batch job 1667690

44)

```
nano star_46_1_pcomp.sh

#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/46_1.fastq.trim.fq _TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/46_1.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pcomp/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pcomp/Pcomp.GFFannotation.fixed_transcript.gff \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pcomp/IndividualStarRuns/46_1.fastq.trim.fq 

sbatch star_46_1_pcomp.sh
```
Submitted batch job 1667692


45)

```
nano star_47_2_pcomp.sh

#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/47_2.fastq.trim.fq _TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/47_2.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pcomp/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pcomp/Pcomp.GFFannotation.fixed_transcript.gff \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pcomp/IndividualStarRuns/47_2.fastq.trim.fq 

sbatch star_47_2_pcomp.sh
```
Submitted batch job 1667693


46)

```
nano star_48_2_pcomp.sh

#!/bin/bash
#SBATCH --mem=10G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/48_2.fastq.trim.fq _TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/48_2.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pcomp/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pcomp/Pcomp.GFFannotation.fixed_transcript.gff \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pcomp/IndividualStarRuns/48_2.fastq.trim.fq 

sbatch star_48_2_pcomp.sh
```
Submitted batch job 1667694


Ran all of these scripts, all of them did not work. They don't fail but they run out of memory after first pass. I upped the memory to 100 GB and am currently rerunning 10_2 with upped memory to see if that works

lol also it wasn't linked the the _tmp in any of these scripts 

Attempt 2) all samples 

Going to run all samples again with full memory and see what happens 

```
nano star_align_batch_script_pcomp_test_upMem.sh
```

```
#!/bin/bash
#SBATCH --mem=64GB

module load STAR/2.5.3a-foss-2016b

FILENAME=$(ls -l /data/putnamlab/jillashey/Francois_data/data/trimmed/ | sed -n $((${SLURM_ARRAY_TASK_ID}+1))p | awk '{print $9}' )

echo $FILENAME

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/${FILENAME}_TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/${FILENAME} \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pcomp/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pcomp/Pcomp.GFFannotation.fixed_transcript.gff
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pcomp_test_upMem/${FILENAME}.
```

```
sbatch --array 1-64 star_align_batch_script_pcomp_test_upMem.sh
```

Submitted batch job 1667865

Errors:
EXITING because of fatal ERROR: could not make run-time genome directory directory: ./_STARgenome/
SOLUTION: please check the path and writing permissions 

Aug 06 11:25:30 ...... FATAL ERROR, exiting
/var/spool/slurmd/job1667882/slurm_script: line 18: --twopassMode: command not found

Going to just go with one sample and test how the script with full mem works


```
nano star_align_batch_script_pcomp_test_upMem_48_2.sh
```

```
#!/bin/bash
#SBATCH --mem=64GB

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/48_2_test_TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/48_2.fastq.trim.fq  \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pcomp/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pcomp/Pcomp.GFFannotation.fixed_transcript.gff
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pcomp_test_upMem/48_2_test.
```

```
sbatch star_align_batch_script_pcomp_test_upMem_48_2.sh
```

Submitted batch job 1667961

1667989_1


Going to try to take out the GTF file and see how it goes

```
nano star_align_16_2_test_noGFF.sh
```

```
#!/bin/bash
#SBATCH --mem=64GB

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/16_2_test_TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/16_2.fastq.trim.fq  \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pcomp/ \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pcomp_test_upMem/16_2_test_noGFF.
```

```
sbatch star_align_16_2_test_noGFF.sh
```
Submitted batch job 1669002
Similar issue to mcap, segmentation fault. Going to try putting it all on one line and with sample 13_2

```
nano star_align_13_2_pcomp_test_noGFF.sh
```

```
#!/bin/bash
#SBATCH --mem=64GB

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads --quantMode TranscriptomeSAM --outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/13_2_pcomp_test_TMP --readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/13_2.fastq.trim.fq --genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pcomp/ --twopassMode Basic --twopass1readsN -1 --outStd Log BAM_Unsorted BAM_Quant --outSAMtype BAM Unsorted SortedByCoordinate --outReadsUnmapped Fastx --outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pcomp_test_upMem/13_2_test_noGFF.
```

```
sbatch star_align_13_2_pcomp_test_noGFF.sh
```

Submitted batch job 1669027
Still says segmentation fault...maybe it is something to do with gff file still impacting it through the genome index
