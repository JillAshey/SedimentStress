### STAR pdam genome

Using annotation file (GFFS) and fasta file (scaffold) from http://pdam.reefgenomics.org/download/ 

a) Make genome index


```
module load STAR/2.5.3a-foss-2016b

STAR --runThreadN 10 --runMode genomeGenerate --genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pdam_rgGFF/ --genomeFastaFiles /data/putnamlab/jillashey/genome/Pdam/pdam_scaffolds.fasta --sjdbGTFfile /data/putnamlab/jillashey/genome/Pdam/pdam_annotation.gff3
```

Not working with that gff; keeps giving me the error: terminate called after throwing an instance of 'std::out_of_range'
  what():  vector::_M_range_check: __n (which is 0) >= this->size() (which is 0)
Aborted
So I put transcript ID in the gff file by taking it from the ID in the attribute column. If this doesnt work, might have to cut some things off the end of transcript_id


```
module load STAR/2.5.3a-foss-2016b

STAR --runThreadN 10 --runMode genomeGenerate --genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pdam_rgGFF/ --genomeFastaFiles /data/putnamlab/jillashey/genome/Pdam/pdam_scaffolds.fasta --sjdbGTFfile /data/putnamlab/jillashey/genome/Pdam/pdam_annotation_AddTranscript_id_fixed.gtf
```


b) Align reads 

Testing a few first 

```
nano star_align_script_pdam_rgGFF.sh
```

```
#!/bin/bash
#SBATCH --mem=100G

module load STAR/2.5.3a-foss-2016b

FILENAME=$(ls -l /data/putnamlab/jillashey/Francois_data/data/trimmed/ | sed -n $((${SLURM_ARRAY_TASK_ID}+1))p | awk '{print $9}' )

STAR --runMode alignReads --quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pdam_rgGFF/${FILENAME}_TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/${FILENAME} \ 
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pdam_rgGFF/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript Parent \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pdam/pdam_annotation_AddTranscript_id_fixed.gtf \
--twopassMode Basic \ 
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \ 
--outSAMtype BAM Unsorted SortedByCoordinate \ 
--outReadsUnmapped Fastx \ 
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pdam_rgGFF/${FILENAME}.

sbatch --array 1-3 star_align_script_pdam_rgGFF.sh

```

Submitted batch job 1672307

Alright they ran well. I'm going to submit them all now, but not as an array job. Because the samples are not clearly marked in a 1,2,3 order, they may be just running randomly. Not sure though

```
sbatch star_align_script_pdam_rgGFF.sh
```

Submitted batch job 1672333

okay well apparently that sbatch job didnt work. says it ran, but it did not. Going to submit as array job and see what happens 

```
sbatch --array 1-64 star_align_script_pdam_rgGFF.sh
```

Only one sample left, will run individually

```
nano star_align_27_1_pdam_rgGFF.sh

#!/bin/bash
#SBATCH --mem=64G

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads --quantMode TranscriptomeSAM --outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/27_1_TMP --readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/27_1.fastq.trim.fq --genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_pdam_rgGFF/ --sjdbGTFfeatureExon exon --sjdbGTFtagExonParentTranscript Parent --sjdbGTFfile /data/putnamlab/jillashey/genome/Pdam/pdam_annotation_AddTranscript_id_fixed.gtf --twopassMode Basic --twopass1readsN -1 --outStd Log BAM_Unsorted BAM_Quant --outSAMtype BAM Unsorted SortedByCoordinate --outReadsUnmapped Fastx --outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pdam/27_1_test_noGFF.

sbatch star_align_27_1_pdam_rgGFF.sh
```
Submitted batch job 1672890
