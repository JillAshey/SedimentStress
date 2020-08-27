# M. capitata 

a-c) in other files but will add here later

d) Generate genome index 

Attempt 1

```
module load STAR/2.5.3a-foss-2016b

STAR --runThreadN 10 --runMode genomeGenerate --genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_mcap/ --genomeFastaFiles /data/putnamlab/jillashey/genome/Mcap/Mcap.genome_assembly.fa --sjdbGTFfile /data/putnamlab/jillashey/genome/Mcap/Mcap.GFFannotation.gff
```

Got this error:
Fatal INPUT FILE error, no exon lines in the GTF file: /data/putnamlab/jillashey/genome/Mcap/Mcap.GFFannotation.gff
Solution: check the formatting of the GTF file, it must contain some lines with exon in the 3rd column.
          Make sure the GTF file is unzipped.
          If exons are marked with a different word, use --sjdbGTFfeatureExon .
Aug 04 13:55:12 ...... FATAL ERROR, exiting

Attempt 2

Adding the ```--sjdbGTFfeatureExon .``` argument because no exons are in the id column. Might be labeled as something else...?

```
module load STAR/2.5.3a-foss-2016b

STAR --runThreadN 10 --runMode genomeGenerate --genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_mcap/ --genomeFastaFiles /data/putnamlab/jillashey/genome/Mcap/Mcap.genome_assembly.fa --sjdbGTFfile /data/putnamlab/jillashey/genome/Mcap/Mcap.GFFannotation.gff --sjdbGTFfeatureExon .
```

Same error as before...Error happening because the features (or id) column of Mcap annotation only has CSD, introns, and genes in it. STAR needs exon features to run. 
Possible solution: if you have no splices (i.e. all transcript are single-exon), you do not need to use the annotations (GTF file) at all. Or, indeed, you can use a non-splice-aware aligner like bwa or bowtie.
If you want to count reads per gene, you can rename all features in your GTF as "exon".

Attempt 3

Feature attributes in gff file were gene, intron, CDS; replaced all with exon

```
module load STAR/2.5.3a-foss-2016b

STAR --runThreadN 10 --runMode genomeGenerate --genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_mcap/ --genomeFastaFiles /data/putnamlab/jillashey/genome/Mcap/Mcap.genome_assembly.fa --sjdbGTFfile /data/putnamlab/jillashey/genome/Mcap/Mcap.GFFannotation.gff
```

Worked!

Attempt 4

Need to change features column so that transcript_id="" has no quotes - no this might be okay, not clear if its best to use transcript_id "" or transcript_id=


e) Align reads 

Try with one sample 

```
nano test_staralign_mcap.sh
```

```
#!/bin/bash
#SBATCH --mem=100GB

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/test_TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/12_1.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_mcap/ \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript ID \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Mcap/Mcap.GFFannotation.gff \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_mcap/12_1_test.
```

```
sbatch test_staralign_mcap.sh
```
Submitted batch job 1667898


Having some issues potentially with the Mcap annotation file, so I am going to take out the GFF annotation file and run it without it (like I did with Plutea)

```
nano star_align_mcap_16_2_test_noGFF.sh
```

```
#!/bin/bash
#SBATCH --mem=100GB

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/16_2_test_noGFF_TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/16_2.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_mcap/ \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_mcap/16_2_test_noGFF.
```

```
sbatch star_align_mcap_16_2_test_noGFF.sh
```

Submitted batch job 1669001
Gave me a segmentation error. Going to try putting it all on one line and with sample 13_2

```
nano star_align_mcap_13_2_test_noGFF.sh
```

```
#!/bin/bash
#SBATCH --mem=100GB

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads --quantMode TranscriptomeSAM --outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/13_2_test_noGFF_TMP --readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/13_2.fastq.trim.fq --genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_mcap/ --twopassMode Basic --twopass1readsN -1 --outStd Log BAM_Unsorted BAM_Quant --outSAMtype BAM Unsorted SortedByCoordinate --outReadsUnmapped Fastx --outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_mcap_noGFF/13_2_test_noGFF.
```

```
sbatch star_align_mcap_13_2_test_noGFF.sh
```
Submitted batch job 1669026

Still says segmentation fault...maybe it is something to do with gff file still impacting it through the genome index




Even though I was able to build genome index, I am still having issues with STAR align for mcap. Took a look at the original and 'fixed' (replaced ids with exon only in bash) annotation files, and in changing ids to exon only, it also changed the gene_id in gene column to exon_id. Not sure if that is the problem, but I changed id to exon in R so it wouldnt mess with the gene columns. Now I will build an index using the updated annotation file and then attempt to use it for alignment 

Feature attributes in gff file were gene, intron, CDS; replaced all with exon

```
module load STAR/2.5.3a-foss-2016b

STAR --runThreadN 10 --runMode genomeGenerate --genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_mcap_all_exons/ --genomeFastaFiles /data/putnamlab/jillashey/genome/Mcap/Mcap.genome_assembly.fa --sjdbGTFfile /data/putnamlab/jillashey/genome/Mcap/Mcap.GFFannotation.fixed_all_exons.gff
```

Worked! Now I need to test it on a sample with STAR alignment

```
nano star_align_mcap_13_2_test_noGFF_allexon.sh
```

```
#!/bin/bash
#SBATCH --mem=100GB

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads --quantMode TranscriptomeSAM --outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/13_2_test_noGFF_allexon_TMP --readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/13_2.fastq.trim.fq --genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_mcap_all_exons/ --twopassMode Basic --twopass1readsN -1 --outStd Log BAM_Unsorted BAM_Quant --outSAMtype BAM Unsorted SortedByCoordinate --outReadsUnmapped Fastx --outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_mcap_noGFF_all_exons/13_2_test_noGFF_allexon.
```

```
sbatch star_align_mcap_13_2_test_noGFF_allexon.sh
```
Submitted batch job 1669065
Still says segmentation fault...

Made some more adjustments with the GFF annotation file in R. I took out the " " around transcript_id and gene_id, set gene names as transcript_id= and gene_id=, took extra spaces and characters out, and changed all ids (intron, CDS, gene) to exon so it could be processed by STAR. ow I will build an index using the updated annotation file and then attempt to use it for alignment !!! Again!

```
module load STAR/2.5.3a-foss-2016b

STAR --runThreadN 10 --runMode genomeGenerate --genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_mcap_all_exons_no_spaces/ --genomeFastaFiles /data/putnamlab/jillashey/genome/Mcap/Mcap.genome_assembly.fa --sjdbGTFfile /data/putnamlab/jillashey/genome/Mcap/Mcap.GFFannotation.fixed_all_exons_no_spaces.gff
```

Worked! Now I need to test it on a sample with STAR alignment

```
nano star_align_mcap_13_2_test_noGFF_allexon_nospaces.sh
```

```
#!/bin/bash
#SBATCH --mem=100GB

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads --quantMode TranscriptomeSAM --outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/13_2_test_noGFF_allexon_nospaces_TMP --readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/13_2.fastq.trim.fq --genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_mcap_all_exons_no_spaces/ --twopassMode Basic --twopass1readsN -1 --outStd Log BAM_Unsorted BAM_Quant --outSAMtype BAM Unsorted SortedByCoordinate --outReadsUnmapped Fastx --outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_mcap_noGFF_all_exons/13_2_test_noGFF_allexon_nospaces.
```

```
sbatch star_align_mcap_13_2_test_noGFF_allexon_nospaces.sh
```

Submitted batch job 1669597
Segmentation error.............

Try once more with another sample 

```
nano star_align_mcap_31_2_test_noGFF_allexon_nospaces.sh
```

```
#!/bin/bash
#SBATCH --mem=100GB

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads --quantMode TranscriptomeSAM --outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/31_2_test_noGFF_allexon_nospaces_TMP --readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/31_2.fastq.trim.fq --genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_mcap_all_exons_no_spaces/ --twopassMode Basic --twopass1readsN -1 --outStd Log BAM_Unsorted BAM_Quant --outSAMtype BAM Unsorted SortedByCoordinate --outReadsUnmapped Fastx --outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_mcap_noGFF_all_exons/31_2_test_noGFF_allexon_nospaces.
```

```
sbatch star_align_mcap_31_2_test_noGFF_allexon_nospaces.sh
```
Submitted batch job 1669620

Trying with all samples

```
nano star_align_batch_script_mcap.sh
```

```
#!/bin/bash
#SBATCH --mem=64G

module load STAR/2.5.3a-foss-2016b

FILENAME=$(ls -l /data/putnamlab/jillashey/Francois_data/data/trimmed/ | sed -n $((${SLURM_ARRAY_TASK_ID}+1))p | awk '{print $9}' )

STAR --runMode alignReads --quantMode TranscriptomeSAM --outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/${FILENAME}_TMP --readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/${FILENAME} --genomeDir /data/putnamlab/jillashey/Francois_data/output/STAR/GenomeIndex_mcap_all_exons_no_spaces/ --twopassMode Basic --twopass1readsN -1 --outStd Log BAM_Unsorted BAM_Quant --outSAMtype BAM Unsorted SortedByCoordinate\ --outReadsUnmapped Fastx --outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_mcap_noGFF_all_exons_no_spaces/${FILENAME}.

sbatch star_align_batch_script_mcap.sh
```
Submitted batch job 1671808

