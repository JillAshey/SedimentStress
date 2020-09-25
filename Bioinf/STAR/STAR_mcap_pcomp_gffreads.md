Running STAR with annotations files processed by gffreads 

gffreads - with both Mcap and Pcomp

### Mcap

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

sbatch --array=1-64%10 AlignReads_pcomp_gtf.sh
```

Submitted batch job 1762269
