Running STAR with annotations files processed by gffreads 

gffreads - with both Mcap and Pcomp

### Mcap

#### Align reads - STAR

a) Convert gff to gtf 

```
gffread Mcap.GFFannotation.fixed.gff -T -o mcap.annotation.gtf
```

STAR with new mcap gtf file 

b) Generate genome index 

```
module load STAR/2.5.3a-foss-2016b

STAR --runThreadN 10 \
--runMode genomeGenerate \
--genomeDir /data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/GenomeIndex_mcap_gtf \
--genomeFastaFiles /data/putnamlab/jillashey/genome/Mcap/Mcap.genome_assembly.fa \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Mcap/mcap.annotation.gtf
```
Not working...says illegal instruction, code dumped. I'm going to get off VPN / putnam lab node and try on bluewaves

When on bluewaves, I use this version of STAR: STAR/2.5.3a-foss-2016b

Weird. STAR running on bluewaves, but not on putnam lab node. 
Okay generating genome worked!

c) Align reads 

```
mkdir AlignReads_mcap_gtf
cd AlignReads_mcap_gtf
ln -s /data/putnamlab/jillashey/Francois_data/Hawaii/data/trimmed/*trim.fq .

# Doing a test run first 

nano test_12_1_staralign_mcap.sh

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

STAR --runMode alignReads --quantMode TranscriptomeSAM --outTmpDir /data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/test_TMP --readFilesIn 12_1.fastq.trim.fq --genomeDir /data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/GenomeIndex_mcap_gtf/ --sjdbGTFfeatureExon exon --sjdbGTFtagExonParentTranscript ID --sjdbGTFfile /data/putnamlab/jillashey/genome/Mcap/mcap.annotation.gtf --twopassMode Basic --twopass1readsN -1 --outStd Log BAM_Unsorted BAM_Quant --outSAMtype BAM Unsorted SortedByCoordinate --outReadsUnmapped Fastx --outFileNamePrefix /data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/AlignReads_mcap_gtf/12_1_test.

```
Submitted batch job 1761987
Okay that didn't work. I'm going to cut the code down so it doesn't have the sjdb stuff in it, similar to what I did with FL samples. 

```
# Doing a test run first 

nano test_12_1_staralign_mcap.sh

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

STAR --runMode alignReads --quantMode TranscriptomeSAM --outTmpDir /data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/test_TMP --readFilesIn 12_1.fastq.trim.fq --genomeDir /data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/GenomeIndex_mcap_gtf/ --twopassMode Basic --twopass1readsN -1 --outStd Log BAM_Unsorted BAM_Quant --outSAMtype BAM Unsorted SortedByCoordinate --outReadsUnmapped Fastx --outFileNamePrefix 12_1_test.

sbatch test_12_1_staralign_mcap.sh
```
Submitted batch job 1761988
Worked!! Can now run all samples!

```
nano AlignReads_mcap_gtf.sh

#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="align_mcap_out_error"
#SBATCH --output="align_mcap_out_error"

module load STAR/2.5.3a-foss-2016b

F=/data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/AlignReads_mcap_gtf

array1=($(ls $F/*trim.fq))

for i in ${array1[@]}
do
STAR --runMode alignReads --quantMode TranscriptomeSAM --outTmpDir ${i}_TMP --readFilesIn ${i} --genomeDir /data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/GenomeIndex_mcap_gtf/ --twopassMode Basic --twopass1readsN -1 --outStd Log BAM_Unsorted BAM_Quant --outSAMtype BAM Unsorted SortedByCoordinate --outReadsUnmapped Fastx --outFileNamePrefix ${i}.
done

sbatch --array=1-64%10 AlignReads_mcap_gtf.sh
```
Submitted batch job 1762259

For some reason, 13_1 and 48_2 didn't run. going to try running individually 

```
## 13_1
nano test_13_1_staralign_mcap.sh

#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="13_test_align_mcap_out_error"
#SBATCH --output="13_test_align_mcap_out_error"

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads --quantMode TranscriptomeSAM --outTmpDir /data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/test13_TMP --readFilesIn 13_1.fastq.trim.fq --genomeDir /data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/GenomeIndex_mcap_gtf/ --twopassMode Basic --twopass1readsN -1 --outStd Log BAM_Unsorted BAM_Quant --outSAMtype BAM Unsorted SortedByCoordinate --outReadsUnmapped Fastx --outFileNamePrefix 13_1_test.

sbatch test_13_1_staralign_mcap.sh

# Submitted batch job 1762818

## 48_2
nano test_48_2_staralign_mcap.sh

#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="48_2_test_align_mcap_out_error"
#SBATCH --output="48_2_test_align_mcap_out_error"

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads --quantMode TranscriptomeSAM --outTmpDir /data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/test48_2_TMP --readFilesIn 48_2.fastq.trim.fq --genomeDir /data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/GenomeIndex_mcap_gtf/ --twopassMode Basic --twopass1readsN -1 --outStd Log BAM_Unsorted BAM_Quant --outSAMtype BAM Unsorted SortedByCoordinate --outReadsUnmapped Fastx --outFileNamePrefix 48_2_test.

sbatch test_48_2_staralign_mcap.sh

## Submitted batch job 1762819
```


#### Count genes - stringTie


a) Move BAM files to stringTie folder 

```
cd /data/putnamlab/jillashey/Francois_data/Florida/output/STAR/AlignReads_mcap_gtf
mv *Aligned.sortedByCoord.out.bam ../../../stringTie_star/mcap_gtf/BAM
```

b) Assemble and estimate reads 

```
nano stringTie_mcap_gtf_assemble.sh

#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="mcap_assemble_out_error"
#SBATCH --output="mcap_assemble_out"

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/Hawaii/stringTie_star/mcap_gtf/BAM

array1=($(ls $F/*bam))
for i in ${array1[@]}; do
	stringtie -G /data/putnamlab/jillashey/genome/Mcap/mcap.annotation.gtf -e -o ${i}.gtf ${i}
	echo "${i}"
done

sbatch stringTie_mcap_gtf_assemble.sh
```

Submitted batch job 1762814

c) Merge stringTie gtf results 

Usually, I do stringTie with all the samples with each annotation file. Now, I think I'm going to move forward only with samples that are assigned to each gff. So only mcap samples from now on


```
mv 5_1.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf, 5_2.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf, 10_2.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf, 13_1.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf, 13_2.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf, 14_1.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf, 14_2.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf, 15_2.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf, 16_2.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf, 19_2.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf, 20_2.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf, 40_2.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf, 45_1.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf, 45_2.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf, 48_2.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf ../GTF

ls *gtf > mcap_mergelist.txt
cat mcap_mergelist.txt

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

stringtie --merge -p 8 -G /data/putnamlab/jillashey/genome/Mcap/mcap.annotation.gtf -o stringtie_mcap_merged.gtf mcap_mergelist.txt
```
d) Assess assembly quality

```
module load gffcompare/0.11.5-foss-2018b

gffcompare -r /data/putnamlab/jillashey/genome/Mcap/mcap.annotation.gtf -o mcap.merged stringtie_mcap_merged.gtf

 63227 reference transcripts loaded.
  63227 query transfrags loaded.
```

e) Re-estimate assembly 

```
nano stringTie_mcap_re-assemble.sh

#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="re-assemble_mcap_out_error"
#SBATCH --output="re-assemble_mcap_out"

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/Hawaii/stringTie_star/mcap_gtf/BAM


array1=($(ls $F/*bam))
for i in ${array1[@]}
do
stringtie -e -G /data/putnamlab/jillashey/genome/Mcap/mcap.annotation.gtf -o ${i}.merge.gtf ${i}
echo "${i}"
done

sbatch stringTie_mcap_re-assemble.sh
```
Submitted batch job 1762821

```
mv 5_1.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf, 5_2.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf,
10_2.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf,
13_1.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf,
13_2.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf, 14_1.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf, 14_2.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf, 15_2.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf, 16_2.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf, 19_2.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf, 20_2.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf, 40_2.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf, 45_1.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf, 
45_2.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf, 
48_2.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf ../GTF_merge
```

f) Create gene matrix

```
module load StringTie/2.1.1-GCCcore-7.3.0
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/Hawaii/stringTie_star/mcap_gtf/GTF_merge/

array2=($(ls *merge.gtf))

for i in ${array2[@]}
do
echo "${i} $F${i}" >> sample_list_mcap.txt
done

python prepDE.py -g gene_count_mcap_gtf_matrix.csv -i sample_list_mcap.txt
```

g) Secure-copy gene counts onto local computer

```
scp jillashey@bluewaves.uri.edu:/data/putnamlab/jillashey/Francois_data/Hawaii/stringTie_star/mcap_gtf/GTF_merge/gene_count_mcap_gtf_matrix.csv /Users/jillashey/Desktop/Putnamlab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/
```


####mcap - When I tried to run STAR on putnam lab node 
Running on putnam lab node

```
mkdir AlignReads_mcap_gtf
cd AlignReads_mcap_gtf
ln -s /data/putnamlab/jillashey/Francois_data/Hawaii/data/trimmed/*trim.fq .

nano test_12_1_staralign_mcap.sh

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

module load STAR/2.6.1c-foss-2018b

STAR --runMode alignReads --quantMode TranscriptomeSAM --outTmpDir /data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/test_TMP --readFilesIn 12_1.fastq.trim.fq --genomeDir /data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/GenomeIndex_mcap/ --sjdbGTFfeatureExon exon --sjdbGTFtagExonParentTranscript ID --sjdbGTFfile /data/putnamlab/jillashey/genome/Mcap/mcap.annotation.gtf --twopassMode Basic --twopass1readsN -1 --outStd Log BAM_Unsorted BAM_Quant --outSAMtype BAM Unsorted SortedByCoordinate --outReadsUnmapped Fastx --outFileNamePrefix /data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/staralign_mcap/12_1_test.

sbatch test_12_1_staralign_mcap.sh
```

Submitted batch job 10575

Getting an error: Sep 24 16:15:46 ..... started STAR run
Sep 24 16:15:48 ..... loading genome
Sep 24 16:16:44 ..... processing annotations GTF

Fatal INPUT FILE error, no exon lines in the GTF file: /data/putnamlab/jillashey/genome/Mcap/mcap.annotation.gtf
Solution: check the formatting of the GTF file, it must contain some lines with exon in the 3rd column.
          Make sure the GTF file is unzipped.
          If exons are marked with a different word, use --sjdbGTFfeatureExon .

Sep 24 16:16:45 ...... FATAL ERROR, exiting

Maybe take out the gtf 

Not working on Putnam lab node. not sure why? Switching over to bluewaves 



### Pcomp

#### Align reads - STAR

a) Convert gff to gtf 

```
gffread Pcomp.GFFannotation.fixed_transcript.gff -T -o pcomp.annotation.gtf

```

STAR with new mcap gtf file 

b) Generate genome index 

```
module load STAR/2.5.3a-foss-2016b

STAR --runThreadN 10 \
--runMode genomeGenerate \
--genomeDir /data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/GenomeIndex_pcomp_gtf \
--genomeFastaFiles /data/putnamlab/jillashey/genome/Pcomp/Porites_compressa_contigs.fasta \
--sjdbGTFfile /data/putnamlab/jillashey/genome/Pcomp/pcomp.annotation.gtf
```
c) Align reads 

```
mkdir AlignReads_pcomp_gtf
cd AlignReads_pcomp_gtf
ln -s /data/putnamlab/jillashey/Francois_data/Hawaii/data/trimmed/*trim.fq .

# Doing a test run first 

nano test_12_1_staralign_pcomp.sh

#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="test_align_pcomp_out_error"
#SBATCH --output="test_align_pcomp_out_error"

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads --quantMode TranscriptomeSAM --outTmpDir /data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/test_pcomp_TMP --readFilesIn 12_1.fastq.trim.fq --genomeDir /data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/GenomeIndex_pcomp_gtf/ --twopassMode Basic --twopass1readsN -1 --outStd Log BAM_Unsorted BAM_Quant --outSAMtype BAM Unsorted SortedByCoordinate --outReadsUnmapped Fastx --outFileNamePrefix 12_1_test.

sbatch test_12_1_staralign_pcomp.sh

```

Submitted batch job 1762194

```
nano AlignReads_pcomp_gtf.sh

#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="align_pcomp_out_error"
#SBATCH --output="align_pcomp_out_error"

module load STAR/2.5.3a-foss-2016b

F=/data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/AlignReads_pcomp_gtf

array1=($(ls $F/*trim.fq))

for i in ${array1[@]}
do
STAR --runMode alignReads --quantMode TranscriptomeSAM --outTmpDir ${i}_TMP --readFilesIn ${i} --genomeDir /data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/GenomeIndex_pcomp_gtf/ --twopassMode Basic --twopass1readsN -1 --outStd Log BAM_Unsorted BAM_Quant --outSAMtype BAM Unsorted SortedByCoordinate --outReadsUnmapped Fastx --outFileNamePrefix ${i}.
done

sbatch AlignReads_pcomp_gtf.sh
```

Submitted batch job 1762806

#### Count genes - stringTie


a) Move BAM files to stringTie folder 

```
cd /data/putnamlab/jillashey/Francois_data/Florida/output/STAR/AlignReads_pcomp_gtf
mv *Aligned.sortedByCoord.out.bam ../../../stringTie_star/pcomp_gtf/BAM
```

b) Assemble and estimate reads 

```
nano stringTie_pcomp_gtf_assemble.sh

#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="pcomp_assemble_out_error"
#SBATCH --output="pcomp_assemble_out"

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/Hawaii/stringTie_star/pcomp_gtf/BAM

array1=($(ls $F/*bam))
for i in ${array1[@]}; do
	stringtie -G /data/putnamlab/jillashey/genome/Pcomp/pcomp.annotation.gtf -e -o ${i}.gtf ${i}
	echo "${i}"
done

sbatch stringTie_pcomp_gtf_assemble.sh
```
Submitted batch job 1764266

c) Merge stringTie gtf results 

Usually, I do stringTie with all the samples with each annotation file. Now, I think I'm going to move forward only with samples that are assigned to each gff. So only pcomp samples from now on


```
mv 3_1.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf 3_2.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf 12_1.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf 17_1.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf
18_1.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf 24_1.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf 30_1.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf 32_1.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf 32_2.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf 33_1.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf 37_1.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf 37_2.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf 43_1.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf 43_2.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf 44_1.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf 46_1.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf ../GTF

ls *gtf > pcomp_mergelist.txt
cat pcomp_mergelist.txt

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

stringtie --merge -p 8 -G /data/putnamlab/jillashey/genome/Pcomp/pcomp.annotation.gtf -o stringtie_pcomp_merged.gtf pcomp_mergelist.txt

Gives me lots of errors: 
Error: discarding overlapping duplicate transcript feature (16707-17244) with ID=Pcom.g74714
Error: discarding overlapping duplicate transcript feature (90130-90565) with ID=Pcom.g72107
Error: discarding overlapping duplicate transcript feature (67208-68911) with ID=Pcom.g72937
Segmentation fault

```

d) Assess assembly quality

```
module load gffcompare/0.11.5-foss-2018b

gffcompare -r /data/putnamlab/jillashey/genome/Pcomp/pcomp.annotation.gtf -o pcomp.merged stringtie_pcomp_merged.gtf

Still more errors
Error: discarding overlapping duplicate transcript feature (5380-6705) with ID=Pcom.g53954
Error: discarding overlapping duplicate transcript feature (6706-6912) with ID=Pcom.g53954
Error: discarding overlapping duplicate transcript feature (6913-7589) with ID=Pcom.g53954
  116881 query transfrags loaded.
  
  Not great sensitivity here...
```

e) Re-estimate assembly 

```
nano stringTie_pcomp_re-assemble.sh

#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="re-assemble_pcomp_out_error"
#SBATCH --output="re-assemble_pcomp_out"

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/Hawaii/stringTie_star/pcomp_gtf/BAM


array1=($(ls $F/*bam))
for i in ${array1[@]}
do
stringtie -e -G /data/putnamlab/jillashey/genome/Pcomp/pcomp.annotation.gtf -o ${i}.merge.gtf ${i}
echo "${i}"
done

sbatch stringTie_pcomp_re-assemble.sh
```
Submitted batch job 1764269

```
mv 3_1.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf, 3_2.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf, 12_1.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf, 17_1.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf,
18_1.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf, 24_1.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf, 30_1.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf, 32_1.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf, 32_2.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf, 33_1.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf, 37_1.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf, 37_2.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf, 43_1.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf, 43_2.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf, 44_1.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf, 46_1.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf ../GTF_merge

```

f) Create gene matrix

```
module load StringTie/2.1.1-GCCcore-7.3.0
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/Hawaii/stringTie_star/pcomp_gtf/GTF_merge/

array2=($(ls *merge.gtf))

for i in ${array2[@]}
do
echo "${i} $F${i}" >> sample_list_pcomp.txt
done

python prepDE.py -g gene_count_pcomp_gtf_matrix.csv -i sample_list_pcomp.txt
```

g) Secure-copy gene counts onto local computer

```
scp jillashey@bluewaves.uri.edu:/data/putnamlab/jillashey/Francois_data/Hawaii/stringTie_star/pcomp_gtf/GTF_merge/gene_count_pcomp_gtf_matrix.csv /Users/jillashey/Desktop/Putnamlab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/
```