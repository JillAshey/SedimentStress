## Pipeline for FL species - sediment stress 

Species: _Acropora cervicornis_, _Montestraea cavernosa_, and _Orbicella faveolata_

Note: some of the paths may not be correct, as I have made new directories and moved files around since running this pipeline.

### 1) Check file integrity 

a) Count all files to make sure all downloaded

```
ls -1 | wc -l
```

b) Verify data integrity with md5sum

```
sbatch md5sum *.txt.gz > checkmd5.md5
md5sum -c checkmd5.md5
```
Should output 'OK' next to each file name 

c) Count number of reads per file 

Some files have different @___. There are: HISEQ, HWI

```
sbatch zgrep -c "HISEQ" *txt.gz
sbatch zgrep -c "HWI" *txt.gz
```
Submitted batch job 1689855 - for HISEQ2 and 1689856 - for HWI

```
# @HISEQ read counts
18_T33_Of_VLL.txt.gz:7,191,417
19_T33_Ac_WK.txt.gz:19,665,420
20_T12_Mc_PWC.txt.gz:18,018,163
21_T33_Mc_EOU.txt.gz:15,475,674
24_T12_Ac_FM.txt.gz:13,786,799
26_T12_Of_WCL.txt.gz:13,436,090
29_T23_Mc_PND.txt.gz:0
30_T23_Of_RPG.txt.gz:0
31_T22_Ac_UV.txt.gz:0
32_T22_Of_EVR.txt.gz:0
34_T22_Mc_SVS.txt.gz:12,609,087
35_T43_Ac_MT.txt.gz:23,846,590
36_T43_Of_JJN.txt.gz:13,551,094
37_T13_Ac_ML.txt.gz:14,171,723
38_T23_Ac_IN.txt.gz:21,348,241
39_T13_Mc_FJE.txt.gz:16,994,083
40_T13_Of_GWS.txt.gz:14,510,120
gzip: 43_ctl3_Of_JVP_2.txt.gz: unexpected end of file
43_ctl3_Of_JVP_2.txt.gz:11,461,086
47_T31_Ac_JB.txt.gz:14,229,140
48_T31_Of_JNO.txt.gz:12,280,565
49_T31_Mc_SWQ.txt.gz:14,934,691
50_T21_Of_YZB.txt.gz:0
51_T42_Of_UOF.txt.gz:0
52_T11_Ac_II.txt.gz:0
53_T21_Ac_NH.txt.gz:0
54_T42_Ac_JQ.txt.gz:0
55_T32_Mc_TWP.txt.gz:0
56_T42_Mc_JAW.txt.gz:0
57_T32_Ac_NM.txt.gz:0
58_T21_Mc_EAH.txt.gz:0
59_T11_Of_TQP.txt.gz:0
60_T32_Of_WXY.txt.gz:0
61_T11_Mc_RAP.txt.gz:0

# @HWI read counts
18_T33_Of_VLL.txt.gz:0
19_T33_Ac_WK.txt.gz:0
20_T12_Mc_PWC.txt.gz:0
21_T33_Mc_EOU.txt.gz:0
24_T12_Ac_FM.txt.gz:0
26_T12_Of_WCL.txt.gz:0
29_T23_Mc_PND.txt.gz:15,499,497
30_T23_Of_RPG.txt.gz:8,218,881
31_T22_Ac_UV.txt.gz:13,317,263
32_T22_Of_EVR.txt.gz:12,667,861
33_T43_Mc_RFV.fastq:18,085,652
34_T22_Mc_SVS.txt.gz:0
35_T43_Ac_MT.txt.gz:0
36_T43_Of_JJN.txt.gz:0
37_T13_Ac_ML.txt.gz:0
38_T23_Ac_IN.txt.gz:0
39_T13_Mc_FJE.txt.gz:0
40_T13_Of_GWS.txt.gz:0
gzip: 43_ctl3_Of_JVP_2.txt.gz: unexpected end of file
43_ctl3_Of_JVP_2.txt.gz:0
47_T31_Ac_JB.txt.gz:0
48_T31_Of_JNO.txt.gz:0
49_T31_Mc_SWQ.txt.gz:0
50_T21_Of_YZB.txt.gz:15,084,847
51_T42_Of_UOF.txt.gz:14,442,222
52_T11_Ac_II.txt.gz:17,303,453
53_T21_Ac_NH.txt.gz:10,306,635
54_T42_Ac_JQ.txt.gz:23,542,186
55_T32_Mc_TWP.txt.gz:19,313,398
56_T42_Mc_JAW.txt.gz:11,743,795
57_T32_Ac_NM.txt.gz:11,679,142
58_T21_Mc_EAH.txt.gz:20,792,838
59_T11_Of_TQP.txt.gz:19,297,200
60_T32_Of_WXY.txt.gz:19,349,770
61_T11_Mc_RAP.txt.gz:16,873,902
```

### 2) Run FastQC

a) Make folders for raw FastQC results and scripts
```
cd Francois_data/Florida
mkdir fastqc_results/raw
mkdir scripts
cd scripts
```

b) Write script for checking quality with FastQC and submit as job on bluewaves

The fastqc code was not working with files as .txt.gz, so changed them to .fastq.gz files

```
nano fastqc_raw.sh

#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Francois_data/Florida/scripts
#SBATCH --error="fastqc_out_raw_error"
#SBATCH --output="fastqc_out_raw"

module load FastQC/0.11.8-Java-1.8 

for file in /data/putnamlab/jillashey/Francois_data/Florida/data/raw/*fastq.gz
do
fastqc $file --outdir /data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/raw
done

sbatch fastqc_raw.sh
```
Submitted batch job 1689863

c) Make sure all files were processed

```
cd fastqc_results/raw
ls -1 | wc -l 
```

### 3) Run MultiQC

a) Make folders for raw MultiWC results

```
cd Francois_data/Florida
mkdir multiqc_results/raw
```

b) Run MultiQC. Pretty fast, so don't need to submit job for it 

```
module load MultiQC/1.7-foss-2018b-Python-2.7.15
multiqc /data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/raw/*fastqc.zip -o /data/putnamlab/jillashey/Francois_data/Florida/multiqc_results/raw
```

![](https://raw.githubusercontent.com/JillAshey/SedimentStress/master/Images/QC/raw/FL/fastqc_sequence_counts_plot.png?token=APHKO34QJ2M3FJUIXHSFO3C7J72OU)


![](https://raw.githubusercontent.com/JillAshey/SedimentStress/master/Images/QC/raw/FL/fastqc_per_sequence_gc_content_plot.png?token=APHKO32U5Q4QGJDBDF4K2OS7J72TM)

![](https://raw.githubusercontent.com/JillAshey/SedimentStress/master/Images/QC/raw/FL/fastqc_per_base_sequence_quality_plot.png?token=APHKO337A3DIJYZWAJWWJK27J72U4)

![](https://raw.githubusercontent.com/JillAshey/SedimentStress/master/Images/QC/raw/FL/fastqc_overrepresented_sequencesi_plot.png?token=APHKO36757WLQVUOXQ5AHV27J72WE)

![](https://raw.githubusercontent.com/JillAshey/SedimentStress/master/Images/QC/raw/FL/fastqc_adapter_content_plot.png?token=APHKO3352RTTC2C3GQBEU3C7J72XG)

c) Copy MultiQC files to local computer

```
scp -r jillashey@bluewaves.uri.edu:/data/putnamlab/jillashey/Francois_data/Florida/multiqc_results/raw/multiqc_data /Users/jillashey/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/QC/raw

scp jillashey@bluewaves.uri.edu:/data/putnamlab/jillashey/Francois_data/Florida/multiqc_results/raw/multiqc_report.html /Users/jillashey/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/QC/raw
```

### 4) Trim reads with Trimmomatic 

a) Make trimmed reads folder in all other results folders 

```
cd Francois_data/Florida
mkdir data/trimmed fastqc_results/trimmed multiqc_results/trimmed
```

b) Unzip fastqc files 

Trimmomatic can't process zipped files 

```
cd data/raw
sbatch gunzip *fastq.gz
```

Unable to gunzip 43_ctl3_Of_JVP_2 because the file ends abruptly in the middle of a line

```
gzip: 43_ctl3_Of_JVP_2.fastq.gz: unexpected end of file
```

c) Write script for Trimmomatic and run on bluewaves

```
cd scripts
nano trimmomatic.sh

#!/bin/bash
#SBATCH -t 48:00:00
#SBATCH --nodes=1 --ntasks-per-node=2
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Francois_data/Florida/scripts
#SBATCH --error="trimmomatic_out_error"
#SBATCH --output="trimmomatic_out"

module load Trimmomatic/0.38-Java-1.8

for file in /data/putnamlab/jillashey/Francois_data/Florida/data/raw/*fastq
do
        java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.38.jar SE -phred33 $file $file.trim.fq ILLUMINACLIP:/data/putnamlab/jillashey/Francois_data/Florida/data/Illumina_adapter_reads_PE_SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20 >> TrimmedAmount.txt
done

sbatch trimmomatic.sh
```
Submitted batch job 1690125

d) Move trimmed files to their own folder 

```
mv *trim.fq ../trimmed
```

### 5) Check quality of trimmed files 

a) Check number of files 

```
ls -1 | wc -l
```

b) Check number of reads

```
sbatch zgrep -c "HISEQ" *trim.fq
sbatch zgrep -c "HWI" *trim.fq

@HISEQ
18_T33_Of_VLL.fastq.trim.fq:6,564,052
19_T33_Ac_WK.fastq.trim.fq:19,493,135
20_T12_Mc_PWC.fastq.trim.fq:17,824,242
21_T33_Mc_EOU.fastq.trim.fq:15,249,255
24_T12_Ac_FM.fastq.trim.fq:13,569,940
26_T12_Of_WCL.fastq.trim.fq:13,252,661
29_T23_Mc_PND.fastq.trim.fq:0
30_T23_Of_RPG.fastq.trim.fq:0
31_T22_Ac_UV.fastq.trim.fq:0
32_T22_Of_EVR.fastq.trim.fq:0
33_T43_Mc_RFV.fastq.trim.fq:0
34_T22_Mc_SVS.fastq.trim.fq:12,419,796
35_T43_Ac_MT.fastq.trim.fq:23,624,545
36_T43_Of_JJN.fastq.trim.fq:13,119,260
37_T13_Ac_ML.fastq.trim.fq:14,089,255
38_T23_Ac_IN.fastq.trim.fq:21,204,725
39_T13_Mc_FJE.fastq.trim.fq:16,844,137
40_T13_Of_GWS.fastq.trim.fq:14,369,470
47_T31_Ac_JB.fastq.trim.fq:14,014,406
48_T31_Of_JNO.fastq.trim.fq:12,042,565
49_T31_Mc_SWQ.fastq.trim.fq:14,519,809
50_T21_Of_YZB.fastq.trim.fq:0
51_T42_Of_UOF.fastq.trim.fq:0
52_T11_Ac_II.fastq.trim.fq:0
53_T21_Ac_NH.fastq.trim.fq:0
54_T42_Ac_JQ.fastq.trim.fq:0
55_T32_Mc_TWP.fastq.trim.fq:0
56_T42_Mc_JAW.fastq.trim.fq:0
57_T32_Ac_NM.fastq.trim.fq:0
58_T21_Mc_EAH.fastq.trim.fq:0
59_T11_Of_TQP.fastq.trim.fq:0
60_T32_Of_WXY.fastq.trim.fq:0
61_T11_Mc_RAP.fastq.trim.fq:0

@HWI
18_T33_Of_VLL.fastq.trim.fq:0
19_T33_Ac_WK.fastq.trim.fq:0
20_T12_Mc_PWC.fastq.trim.fq:0
21_T33_Mc_EOU.fastq.trim.fq:0
24_T12_Ac_FM.fastq.trim.fq:0
26_T12_Of_WCL.fastq.trim.fq:0
29_T23_Mc_PND.fastq.trim.fq:15,429,503
30_T23_Of_RPG.fastq.trim.fq:8,205,037
31_T22_Ac_UV.fastq.trim.fq:13,262,941
32_T22_Of_EVR.fastq.trim.fq:12,482,192
33_T43_Mc_RFV.fastq.trim.fq:17,860,108
34_T22_Mc_SVS.fastq.trim.fq:0
35_T43_Ac_MT.fastq.trim.fq:0
36_T43_Of_JJN.fastq.trim.fq:0
37_T13_Ac_ML.fastq.trim.fq:0
38_T23_Ac_IN.fastq.trim.fq:0
39_T13_Mc_FJE.fastq.trim.fq:0
40_T13_Of_GWS.fastq.trim.fq:0
47_T31_Ac_JB.fastq.trim.fq:0
48_T31_Of_JNO.fastq.trim.fq:0
49_T31_Mc_SWQ.fastq.trim.fq:0
50_T21_Of_YZB.fastq.trim.fq:14,983,598
51_T42_Of_UOF.fastq.trim.fq:14,396,383
52_T11_Ac_II.fastq.trim.fq:17,145,972
53_T21_Ac_NH.fastq.trim.fq:10,130,277
54_T42_Ac_JQ.fastq.trim.fq:23,414,821
55_T32_Mc_TWP.fastq.trim.fq:19,235,526
56_T42_Mc_JAW.fastq.trim.fq:11,663,729
57_T32_Ac_NM.fastq.trim.fq:11,628,199
58_T21_Mc_EAH.fastq.trim.fq:20,509,453
59_T11_Of_TQP.fastq.trim.fq:19,238,856
60_T32_Of_WXY.fastq.trim.fq:19,279,940
61_T11_Mc_RAP.fastq.trim.fq:16,801,355
```
Submitted batch job 1690400 for HISEQ and 1690401 for HWI

c) Run FastQC on trimmed data

```
nano fastqc_trimmed.sh

#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Francois_data/Florida/scripts
#SBATCH --error="fastqc_out_trimmed_error"
#SBATCH --output="fastqc_out_trimmed"

module load FastQC/0.11.8-Java-1.8 

for file in /data/putnamlab/jillashey/Francois_data/Florida/data/trimmed/*trim.fq
do
fastqc $file --outdir /data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/trimmed
done

sbatch fastqc_trimmed.sh
```
Submitted batch job 1690402

d) Run MultiQC on trimmed data and download files to local computer

```
module load MultiQC/1.7-foss-2018b-Python-2.7.15
multiqc /data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/trimmed/*fastqc.zip -o /data/putnamlab/jillashey/Francois_data/Florida/multiqc_results/trimmed

scp -r jillashey@bluewaves.uri.edu:/data/putnamlab/jillashey/Francois_data/Florida/multiqc_results/trimmed/multiqc_data /Users/jillashey/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/QC/trimmed

scp jillashey@bluewaves.uri.edu:/data/putnamlab/jillashey/Francois_data/Florida/multiqc_results/trimmed/multiqc_report.html /Users/jillashey/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/QC/trimmed
```

![](https://raw.githubusercontent.com/JillAshey/SedimentStress/master/Images/QC/trimmed/fastqc_sequence_counts_plot.png?token=APHKO32IMCD4U5AQOXNGYMS7J726W)

![](https://raw.githubusercontent.com/JillAshey/SedimentStress/master/Images/QC/trimmed/fastqc_per_sequence_gc_content_plot.png?token=APHKO37KY2NCDO56RVLUQPC7J73CC)

![](https://raw.githubusercontent.com/JillAshey/SedimentStress/master/Images/QC/trimmed/fastqc_per_base_sequence_quality_plot.png?token=APHKO32JTY6U2JJXCY2J3U27J73JW)

![](https://raw.githubusercontent.com/JillAshey/SedimentStress/master/Images/QC/trimmed/fastqc_overrepresented_sequencesi_plot.png?token=APHKO33VGTBNNZAGIKCKSD27J73OY)

![](https://raw.githubusercontent.com/JillAshey/SedimentStress/master/Images/QC/trimmed/fastqc_adapter_content_plot.png?token=APHKO334NFOXWGH7KFVQZRK7J73PY)

Adapter content looks kinda weird...maybe because Francois may have already trimmed some of the samples? But not sure about that yet

### 6) Align reads with STAR

Before running STAR, I added the identifier 'transcript_id=' to the last column of the gff file in R. STAR needs this identifier to run and most of the gffs I used don't have it. Code to edit gff files herexxxxxxx

#### O. fav

a) Generate genome index

```
module load STAR/2.5.3a-foss-2016b

STAR --runThreadN 10 --runMode genomeGenerate --genomeDir /data/putnamlab/jillashey/Francois_data/Florida/output/STAR/GenomeIndex_Ofav --genomeFastaFiles /data/putnamlab/jillashey/genome/Ofav/GCF_002042975.1_ofav_dov_v1_genomic.fna --sjdbGTFfile /data/putnamlab/jillashey/genome/Ofav/GCF_002042975.1_ofav_dov_v1_genomic.gff
```

b) Align reads to genome


```
ln -s /data/putnamlab/jillashey/Francois_data/Florida/data/trimmed/*trim.fq .

nano test_align_script.sh

#!/bin/bash
#SBATCH -t 48:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="Align_Ofav_out_error"
#SBATCH --output="Align_Ofav_out"

module load STAR/2.5.3a-foss-2016b

F=/data/putnamlab/jillashey/Francois_data/Florida/output/STAR/AlignReads_Ofav

array1=($(ls $F/*trim.fq))
for i in ${array1[@]}
do
STAR --runMode alignReads --quantMode TranscriptomeSAM --outTmpDir ${i}_TMP --readFilesIn ${i} --genomeDir /data/putnamlab/jillashey/Francois_data/Florida/output/STAR/GenomeIndex_Ofav --twopassMode Basic --twopass1readsN -1 --outStd Log BAM_Unsorted BAM_Quant --outSAMtype BAM Unsorted SortedByCoordinate --outReadsUnmapped Fastx --outFileNamePrefix ${i}.
done 

sbatch test_align_script.sh 
```
Submitted batch job 1691569

#### A. cerv

a) Generate genome index

Added transcript_id to gene col in gff file

```
module load STAR/2.5.3a-foss-2016b

STAR --runThreadN 10 --runMode genomeGenerate --genomeDir /data/putnamlab/jillashey/Francois_data/Florida/output/STAR/GenomeIndex_Acerv --genomeFastaFiles /data/putnamlab/jillashey/genome/Acerv/Acerv_assembly_v1.0_171209.fasta --sjdbGTFfile /data/putnamlab/jillashey/genome/Acerv/Acerv.GFFannotations.fixed_transcript.gff3

```

b) Align reads to genome

```
mkdir AlignReads_Acerv
cd AlignReads_Acerv
ln -s /data/putnamlab/jillashey/Francois_data/Florida/data/trimmed/*trim.fq .

nano AlignReads_acerv.sh

#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="Align_Acerv_out_error"
#SBATCH --output="Align_Acerv_out"

module load STAR/2.5.3a-foss-2016b

F=/data/putnamlab/jillashey/Francois_data/Florida/output/STAR/AlignReads_Acerv

array1=($(ls $F/*trim.fq))
for i in ${array1[@]}
do
STAR --runMode alignReads --quantMode TranscriptomeSAM --outTmpDir ${i}_TMP --readFilesIn ${i} --genomeDir /data/putnamlab/jillashey/Francois_data/Florida/output/STAR/GenomeIndex_Acerv --twopassMode Basic --twopass1readsN -1 --outStd Log BAM_Unsorted BAM_Quant --outSAMtype BAM Unsorted SortedByCoordinate --outReadsUnmapped Fastx --outFileNamePrefix ${i}.
done 

sbatch AlignReads_acerv.sh 
```

Submitted batch job 1691808

#### M. cav

a) Generate genome index


Added transcript_id to gene col in gff

```
module load STAR/2.5.3a-foss-2016b

STAR --runThreadN 10 --runMode genomeGenerate --genomeDir /data/putnamlab/jillashey/Francois_data/Florida/output/STAR/GenomeIndex_Mcav --genomeFastaFiles /data/putnamlab/jillashey/genome/Mcav/Mcavernosa_July2018.fasta --sjdbGTFfile /data/putnamlab/jillashey/genome/Mcav/Mcavernosa_annotation/Mcav.gff.annotations.fixed_transcript.gff3

```

b) Align reads to genome

```
mkdir AlignReads_Mcav
cd AlignReads_Mcav
ln -s /data/putnamlab/jillashey/Francois_data/Florida/data/trimmed/*trim.fq .

nano AlignReads_mcav.sh

#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="Align_Acerv_out_error"
#SBATCH --output="Align_Acerv_out"

module load STAR/2.5.3a-foss-2016b

F=/data/putnamlab/jillashey/Francois_data/Florida/output/STAR/AlignReads_Mcav

array1=($(ls $F/*trim.fq))
for i in ${array1[@]}
do
STAR --runMode alignReads --quantMode TranscriptomeSAM --outTmpDir ${i}_TMP --readFilesIn ${i} --genomeDir /data/putnamlab/jillashey/Francois_data/Florida/output/STAR/GenomeIndex_Mcav --twopassMode Basic --twopass1readsN -1 --outStd Log BAM_Unsorted BAM_Quant --outSAMtype BAM Unsorted SortedByCoordinate --outReadsUnmapped Fastx --outFileNamePrefix ${i}.
done 

sbatch AlignReads_mcav.sh 
```
Submitted batch job 1691816


### 6) Perform gene counts with stringTie

#### O. fav

```
mkdir stringTie 
cd stringTie
mkdir Acerv Mcav Ofav

# in each species folder 
mkdir BAM GTF GTF_merge
```

a) Move BAM files to stringTie folder 

```
cd /data/putnamlab/jillashey/Francois_data/Florida/output/STAR/AlignReads_Ofav
mv *Aligned.sortedByCoord.out.bam ../../../stringTie/Ofav/BAM
```

b) Assemble and estimate reads 

```
cd stringTie/Ofav/BAM

nano stringTie_ofav_assemble.sh

#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="assemble_ofav_out_error"
#SBATCH --output="assemble_ofav_out"

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/Florida/stringTie/Ofav/BAM/

array1=($(ls $F/*bam))
for i in ${array1[@]}; do
	stringtie -G /data/putnamlab/jillashey/genome/Ofav/GCF_002042975.1_ofav_dov_v1_genomic.gff -e -o ${i}.gtf ${i}
	echo "${i}"
done

sbatch stringTie_ofav_assemble.sh
```

Submitted batch job 1692965


c) Merge stringTie gtf results 

```
mv *gtf ../GTFfiles/

ls *gtf > ofav_mergelist.txt
cat ofav_mergelist.txt

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

stringtie --merge -p 8 -G /data/putnamlab/jillashey/genome/Ofav/GCF_002042975.1_ofav_dov_v1_genomic.gff -o stringtie_ofav_merged.gtf ofav_mergelist.txt
```

d) Assess assembly quality

```
module load gffcompare/0.11.5-foss-2018b

gffcompare -r /data/putnamlab/jillashey/genome/Ofav/GCF_002042975.1_ofav_dov_v1_genomic.gff -o Ofav.merged stringtie_ofav_merged.gtf

37786 reference transcripts loaded.
  5 duplicate reference transcripts discarded.
  37781 query transfrags loaded.
```

e) Re-estimate assembly 

```
nano stringTie_ofav_re-assemble.sh

#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="re-assemble_ofav_out_error"
#SBATCH --output="re-assemble_ofav_out"

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/Florida/stringTie/Ofav/BAM/

array1=($(ls $F/*bam))
for i in ${array1[@]}
do
stringtie -e -G /data/putnamlab/jillashey/genome/Ofav/GCF_002042975.1_ofav_dov_v1_genomic.gff -o ${i}.merge.gtf ${i}
echo "${i}"
done

sbatch stringTie_ofav_re-assemble.sh

mv *merge.gtf ../GTF_merge
```
Submitted batch job 1692972

f) Create gene matrix

```
module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/Florida/stringTie/Ofav/GTF_merge/

array2=($(ls *merge.gtf))

for i in ${array2[@]}
do
echo "${i} $F${i}" >> sample_list_ofav.txt
done

python prepDE.py -g gene_count_ofav_matrix.csv -i sample_list_ofav.txt
```

g) Secure-copy gene counts onto local computer

```
scp jillashey@bluewaves.uri.edu:/data/putnamlab/jillashey/Francois_data/Florida/stringTie/Ofav/GTF_merge/gene_count_ofav_matrix.csv /Users/jillashey/Desktop/Putnamlab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/star/
```

#### A. cerv

a) Move BAM files to stringTie folder 

```
cd /data/putnamlab/jillashey/Francois_data/Florida/output/STAR/AlignReads_Acerv
mv *Aligned.sortedByCoord.out.bam ../../../stringTie/Acerv/BAM
```

b) Assemble and estimate reads 

```
cd stringTie/Acerv/BAM

nano stringTie_acerv_assemble.sh

#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="assemble_acerv_out_error"
#SBATCH --output="assemble_acerv_out"

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/Florida/stringTie/Acerv/BAM/

array1=($(ls $F/*bam))
for i in ${array1[@]}; do
	stringtie -G /data/putnamlab/jillashey/genome/Acerv/Acerv_assembly_v1.0.gff3 -e -o ${i}.gtf ${i}
	echo "${i}"
done

sbatch stringTie_acerv_assemble.sh

```
Submitted batch job 1692969

c) Merge stringTie gtf results 

```
mv *gtf ../GTFfiles/

ls *gtf > acerv_mergelist.txt
cat acerv_mergelist.txt

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

stringtie --merge -p 8 -G /data/putnamlab/jillashey/genome/Acerv/Acerv_assembly_v1.0.gff3 -o stringtie_acerv_merged.gtf acerv_mergelist.txt
```

d) Assess assembly quality

```
module load gffcompare/0.11.5-foss-2018b

gffcompare -r /data/putnamlab/jillashey/genome/Acerv/Acerv_assembly_v1.0.gff3 -o Acerv.merged stringtie_acerv_merged.gtf

  48478 reference transcripts loaded.
  48478 query transfrags loaded.
```

e) Re-estimate assembly 

```
nano stringTie_acerv_re-assemble.sh

#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="re-assemble_acerv_out_error"
#SBATCH --output="re-assemble_acerv_out"

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/Florida/stringTie/Acerv/BAM/

array1=($(ls $F/*bam))
for i in ${array1[@]}
do
stringtie -e -G /data/putnamlab/jillashey/genome/Acerv/Acerv_assembly_v1.0.gff3 -o ${i}.merge.gtf ${i}
echo "${i}"
done

sbatch stringTie_acerv_re-assemble.sh

mv *merge.gtf ../GTF_merge
```
Submitted batch job 1692973


f) Create gene matrix

```
module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/Florida/stringTie/Acerv/GTF_merge/

array2=($(ls *merge.gtf))

for i in ${array2[@]}
do
echo "${i} $F${i}" >> sample_list_acerv.txt
done

python prepDE.py -g gene_count_acerv_matrix.csv -i sample_list_acerv.txt
```

g) Secure-copy gene counts onto local computer

```
scp jillashey@bluewaves.uri.edu:/data/putnamlab/jillashey/Francois_data/Florida/stringTie/Acerv/GTF_merge/gene_count_acerv_matrix.csv /Users/jillashey/Desktop/Putnamlab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/star/
```

#### M. cav

a) Move BAM files to stringTie folder 

```
cd /data/putnamlab/jillashey/Francois_data/Florida/output/STAR/AlignReads_Mcav
mv *Aligned.sortedByCoord.out.bam ../../../stringTie/Mcav/BAM
```

b) Assemble and estimate reads 

```
cd stringTie/Mcav/BAM

nano stringTie_mcav_assemble.sh

#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="assemble_mcav_out_error"
#SBATCH --output="assemble_mcav_out"

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/Florida/stringTie/Mcav/BAM/

array1=($(ls $F/*bam))
for i in ${array1[@]}; do
	stringtie -G /data/putnamlab/jillashey/genome/Mcav/Mcavernosa_annotation/Mcav.gff.annotations.fixed_transcript.gff3 -e -o ${i}.gtf ${i}
	echo "${i}"
done

sbatch stringTie_mcav_assemble.sh

```

Submitted batch job 1693008

Submitted batch job 1693022 -- have to rerun bc accidently ran it with the original gff file instead of the one I fixed and that I used in STAR 


c) Merge stringTie gtf results 

```
mv *gtf ../GTFfiles/

ls *gtf > mcav_mergelist.txt
cat mcav_mergelist.txt

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

stringtie --merge -p 8 -G /data/putnamlab/jillashey/genome/Mcav/Mcavernosa_annotation/Mcav.gff.annotations.fixed_transcript.gff3 -o stringtie_mcav_merged.gtf mcav_mergelist.txt
```

When running with original Mcav gff file: for rows that were only mRNA or gene, I got this error: 

```
Error: invalid feature coordinates (end<start!) at line:
Sc0000079	maker	mRNA	734362	719	.	-	.	ID=Mcavernosa02687-RA;Parent=Mcavernosa02687;Name=Mcavernosa02687-RA;Alias=maker-NewChr8-snap-gene-61.20-mRNA-1;_AED=0.01;_QI=0|0.4|0.33|0.83|0.4|0.5|6|2976|639;_eAED=0.01;

```

When running with fixed Mcav gff file: I got weird NA lines and still errors around mRNA and gene

```
Warning: invalid start coordinate at line:
###			NA	NA				;transcript_id=
Error: invalid feature coordinates (end<start!) at line:
Sc0000070	maker	gene	801834	5681	.	-	.	ID=Mcavernosa17287;Name=Mcavernosa17287;Alias=augustus_masked-NewChr7-processed-gene-107.6;
Error: invalid feature coordinates (end<start!) at line:
Sc0000070	maker	mRNA	801834	5681	.	-	.	ID=Mcavernosa17287-RA;Parent=Mcavernosa17287;Name=Mcavernosa17287-RA;Alias=augustus_masked-NewChr7-processed-gene-107.6-mRNA-1;_AED=0.21;_QI=403|0.6|0.66|0.83|0.6|0.5|6|0|480;_eAED=0.21;;transcript_id=Mcavernosa17287-RA
Warning: invalid start coordinate at line:
###			NA	NA				;transcript_id=
Warning: invalid start coordinate at line:
###			NA	NA				;transcript_id=
Warning: invalid start coordinate at line:
###			NA	NA				;transcript_id=
Warning: invalid start coordinate at line:

```

Not sure what this means...Seems like new gtf file only has exons and transcripts. The error about end<start seems like something is up with the stop and start coords. 

Figured out that some of the gene lengths in the Mcav gff file were negative, which doesn't make any sense because how can a gene be a negative length? So I saved the positive gene lengths only into new Mcav gff file, which I'm now going to run STAR and stringTie on to see if I can fix the errors I was getting above.
________________________________________

Initially, I did not have all the samples to analyze, but Francois recently uploaded them to Tufts server. I copied them over to Bluewaves and I'm going to put those samples through the pipeline. 

```
scp jashey@linux.eecs.tufts.edu:/r/corals/Jill/Sediment_FL/*fastq.gz jillashey@bluewaves.uri.edu:/data/putnamlab/jillashey/Francois_data/Florida/data/raw/raw_include
```


### 1) Check file integrity 

a) Count all files to make sure all downloaded

```
ls -1 | wc -l
```

b) Verify data integrity with md5sum

```
md5sum *.fastq.gz > checkmd5.md5
md5sum -c checkmd5.md5
```
Should output 'OK' next to each file name 

c) Count number of reads per file 

Some files have different @___. There are: HISEQ, HWI

```
sbatch zgrep -c "HISEQ" *fastq.gz
sbatch zgrep -c "HWI" *txt.gz
```
Submitted batch job 1717762 - HISEQ 

```
17_ctl2_Of_ZTH_1.fastq.gz:16,862,085
17_ctl2_Of_ZTH_2.fastq.gz:16,862,085
22_ctl2_Mc_TWF_1.fastq.gz:9,492,297
22_ctl2_Mc_TWF_2.fastq.gz:9492297
23_ctl1_Of_CTX_1.fastq.gz:9,890,631
23_ctl1_Of_CTX_2.fastq.gz:9,890,631
25_ctl1_Ac_GF_1.fastq.gz:14,630,151
25_ctl1_Ac_GF_2.fastq.gz:14,630,151
27_ctl2_Ac_YG_1.fastq.gz:15,708,593
27_ctl2_Ac_YG_2.fastq.gz:15,708,593
28_ctl1_Mc_GBM_1.fastq.gz:18,832,994
28_ctl1_Mc_GBM_2.fastq.gz:18,832,994
41_ctl3_Ac_RN_1.fastq.gz:16,449,541
41_ctl3_Ac_RN_2.fastq.gz:16,449,541
42_ctl3_Mc_MGR_1.fastq.gz:16,522,574
42_ctl3_Mc_MGR_2.fastq.gz:16,522,574
43_ctl3_Of_JVP_1.fastq.gz:16,107,263
43_ctl3_Of_JVP_2.fastq.gz:16,107,263
44_T41_Of_PVT_1.fastq.gz:14,566,583
44_T41_Of_PVT_2.fastq.gz:14,566,583
45_T41_Ac_SC_1.fastq.gz:20,233,577
45_T41_Ac_SC_2.fastq.gz:20,233,577
46_T41_Mc_QYH_1.fastq.gz:13,610,746
46_T41_Mc_QYH_2.fastq.gz:13,610,746
```

Interesting...files 1 and 2 have same # of lines. Does this mean they are just the same files?

```
22
@HISEQ-H454:145:CAL3CANXX:3:1101:4861:1998 1:N:0:AGTCAA
@HISEQ-H454:145:CAL3CANXX:3:1101:4861:1998 2:N:0:AGTCAA

44
@HISEQ-H454:145:CAL3CANXX:3:1101:8803:1999 1:N:0:TAGCTT
@HISEQ-H454:145:CAL3CANXX:3:1101:8803:1999 2:N:0:TAGCTT
```

Appears that it is the same file for each sample just duplicated for some reason. Going to proceed with *_1.fastq.gz samples

### 2) Run FastQC

a) Write script for checking quality with FastQC 

```
nano fastqc_raw_include.sh

#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Francois_data/Florida/scripts
#SBATCH --error="fastqc_out_raw_include_error"
#SBATCH --output="fastqc_out_raw_include"

module load FastQC/0.11.8-Java-1.8 

# These samples are the second batch that Francois uploaded to server

for file in /data/putnamlab/jillashey/Francois_data/Florida/data/raw/raw_include/*1.fastq.gz
do
fastqc $file --outdir /data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/raw/raw_include
done

sbatch fastqc_raw_include.sh
```
Submitted batch job 1717764

b) Make sure all files were processed

```
cd fastqc_results/raw/raw_include
ls -1 | wc -l 
```

### 3) Run MultiQC

a) Make folders for raw MultiWC results

```
cd Francois_data/Florida
mkdir multiqc_results/raw/raw_include
```

b) Run MultiQC. Pretty fast, so don't need to submit job for it 

```
module load MultiQC/1.7-foss-2018b-Python-2.7.15
multiqc /data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/raw/raw_include/*fastqc.zip -o /data/putnamlab/jillashey/Francois_data/Florida/multiqc_results/raw/raw_include
```

![](https://raw.githubusercontent.com/JillAshey/SedimentStress/master/Images/QC/raw/FL/fastqc_sequence_counts_plot_raw_include.png?token=APHKO37I5RUTC3HB24IQL6C7MA7X2)

![](https://raw.githubusercontent.com/JillAshey/SedimentStress/master/Images/QC/raw/FL/fastqc_per_sequence_gc_content_plot_raw_include.png?token=APHKO3452UVGI7MZ3TDBQQ27MA7ZW)

![](https://raw.githubusercontent.com/JillAshey/SedimentStress/master/Images/QC/raw/FL/fastqc_per_base_sequence_quality_plot_raw_include.png?token=APHKO36RJRHJPWLL74JI6NK7MA74C)

![](https://raw.githubusercontent.com/JillAshey/SedimentStress/master/Images/QC/raw/FL/fastqc_overrepresented_sequencesi_plot_raw_include.png?token=APHKO33JUOMCTUU2PQKMSVS7MA75U)

![](https://raw.githubusercontent.com/JillAshey/SedimentStress/master/Images/QC/raw/FL/fastqc_adapter_content_plot_raw_include.png?token=APHKO32GKQCOCEJEPO724LS7MA764)

c) Copy files to local computer 

```
scp jillashey@bluewaves.uri.edu:/data/putnamlab/jillashey/Francois_data/Florida/multiqc_results/raw/raw_include/multiqc_report.html /Users/jillashey/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/QC/raw 

scp -r jillashey@bluewaves.uri.edu:/data/putnamlab/jillashey/Francois_data/Florida/multiqc_results/raw/raw_include/multiqc_data/ /Users/jillashey/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/QC/raw
```

### 4) Trim reads with Trimmomatic 

a) Make trimmed reads folder in all other results folders 

```
cd Francois_data/Florida
mkdir data/trimmed/trim_include fastqc_results/trimmed/trim_include multiqc_results/trimmed/trim_include
```

b) Unzip fastqc files 

Trimmomatic can't process zipped files 

```
cd data/raw/raw_include 
sbatch gunzip *1.fastq.gz
```
Submitted batch job 1717770


c) Write script for Trimmomatic and run on bluewaves

```
cd scripts
nano trimmomatic_include.sh

#!/bin/bash
#SBATCH -t 48:00:00
#SBATCH --nodes=1 --ntasks-per-node=2
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Francois_data/Florida/scripts
#SBATCH --error="trimmomatic_include_out_error"
#SBATCH --output="trimmomatic_include_out"

module load Trimmomatic/0.38-Java-1.8

for file in /data/putnamlab/jillashey/Francois_data/Florida/data/raw/raw_include/*1.fastq
do
        java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.38.jar SE -phred33 $file $file.trim.fq ILLUMINACLIP:/data/putnamlab/jillashey/Francois_data/Florida/data/Illumina_adapter_reads_PE_SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20 >> TrimmedAmount_include.txt
done

sbatch trimmomatic.sh
```
Submitted batch job 1717772

d) Move trimmed files to their own folder 

```
mv *trim.fq ../trimmed/trim_include
```

### 5) Check quality of trimmed files 

a) Check number of files 

```
ls -1 | wc -l
```

b) Check number of reads

```
sbatch zgrep -c "HISEQ" *trim.fq

17_ctl2_Of_ZTH_1.fastq.trim.fq:16,609,827
22_ctl2_Mc_TWF_1.fastq.trim.fq:9,117,346
23_ctl1_Of_CTX_1.fastq.trim.fq:9,672,984
25_ctl1_Ac_GF_1.fastq.trim.fq:14,537,136
27_ctl2_Ac_YG_1.fastq.trim.fq:15,618,761
28_ctl1_Mc_GBM_1.fastq.trim.fq:18,735,353
41_ctl3_Ac_RN_1.fastq.trim.fq:16,334,644
42_ctl3_Mc_MGR_1.fastq.trim.fq:16,398,481
43_ctl3_Of_JVP_1.fastq.trim.fq:15,683,539
44_T41_Of_PVT_1.fastq.trim.fq:14,443,689
45_T41_Ac_SC_1.fastq.trim.fq:20,138,121
46_T41_Mc_QYH_1.fastq.trim.fq:13,494,826
```
Submitted batch job 1719275

c) Run FastQC on trimmed data

```
nano fastqc_trim_include.sh

#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Francois_data/Florida/scripts
#SBATCH --error="fastqc_out_trim_include_error"
#SBATCH --output="fastqc_out_trim_include"

module load FastQC/0.11.8-Java-1.8 

for file in /data/putnamlab/jillashey/Francois_data/Florida/data/trimmed/trim_include/*trim.fq
do
fastqc $file --outdir /data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/trimmed/trim_include
done

sbatch fastqc_trim_include.sh
```
Submitted batch job 1717848

d) Run MultiQC on trimmed data and moved files 

```
module load MultiQC/1.7-foss-2018b-Python-2.7.15
multiqc /data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/trimmed/trim_include/*fastqc.zip -o /data/putnamlab/jillashey/Francois_data/Florida/multiqc_results/trimmed/trim_include 

scp -r jillashey@bluewaves.uri.edu:/data/putnamlab/jillashey/Francois_data/Florida/multiqc_results/trimmed/trim_include/multiqc_data /Users/jillashey/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/QC/trimmed/

scp jillashey@bluewaves.uri.edu:/data/putnamlab/jillashey/Francois_data/Florida/multiqc_results/trimmed/trim_include/multiqc_report.html /Users/jillashey/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/QC/trimmed/
```

### 6) Align reads to genome with STAR

####O. fav

a) Generate genome index - already done from above 

```
# This is the genome index I'll be using 
--genomeDir /data/putnamlab/jillashey/Francois_data/Florida/output/STAR/GenomeIndex_Ofav
```

b) Align reads to genome

```
mkdir AlignReads_Ofav_include
cd AlignReads_Ofav_include
ln -s /data/putnamlab/jillashey/Francois_data/Florida/data/trimmed/trim_include/*trim.fq .

nano AlignReads_ofav_include.sh

#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="Align_Ofav_include_out_error"
#SBATCH --output="Align_Ofav_include_out"

module load STAR/2.5.3a-foss-2016b

F=/data/putnamlab/jillashey/Francois_data/Florida/output/STAR/AlignReads_Ofav_include

array1=($(ls $F/*trim.fq))
for i in ${array1[@]}
do
STAR --runMode alignReads --quantMode TranscriptomeSAM --outTmpDir ${i}_TMP --readFilesIn ${i} --genomeDir /data/putnamlab/jillashey/Francois_data/Florida/output/STAR/GenomeIndex_Ofav --twopassMode Basic --twopass1readsN -1 --outStd Log BAM_Unsorted BAM_Quant --outSAMtype BAM Unsorted SortedByCoordinate --outReadsUnmapped Fastx --outFileNamePrefix ${i}.
done 

sbatch AlignReads_ofav_include.sh 
```
Submitted batch job 1717853

####A. cerv

a) Generate genome index - already done from above 

```
# This is the genome index I'll be using 
--genomeDir /data/putnamlab/jillashey/Francois_data/Florida/output/STAR/GenomeIndex_Acerv 
```

b) Align reads to genome

```
mkdir AlignReads_Acerv_include
cd AlignReads_Acerv_include
ln -s /data/putnamlab/jillashey/Francois_data/Florida/data/trimmed/trim_include/*trim.fq .

nano AlignReads_acerv_include.sh

#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="Align_Acerv_include_out_error"
#SBATCH --output="Align_Acerv_include_out"

module load STAR/2.5.3a-foss-2016b

F=/data/putnamlab/jillashey/Francois_data/Florida/output/STAR/AlignReads_Acerv_include

array1=($(ls $F/*trim.fq))
for i in ${array1[@]}
do
STAR --runMode alignReads --quantMode TranscriptomeSAM --outTmpDir ${i}_TMP --readFilesIn ${i} --genomeDir /data/putnamlab/jillashey/Francois_data/Florida/output/STAR/GenomeIndex_Acerv --twopassMode Basic --twopass1readsN -1 --outStd Log BAM_Unsorted BAM_Quant --outSAMtype BAM Unsorted SortedByCoordinate --outReadsUnmapped Fastx --outFileNamePrefix ${i}.
done 

sbatch AlignReads_acerv_include.sh 
```
Submitted batch job 1717851

####M. cav

a) Generate genome index - already done from above 

```
# This is the genome index I'll be using 
--genomeDir /data/putnamlab/jillashey/Francois_data/Florida/output/STAR/GenomeIndex_Mcav
```

b) Align reads to genome

```
mkdir AlignReads_Mcav_include
cd AlignReads_Mcav_include
ln -s /data/putnamlab/jillashey/Francois_data/Florida/data/trimmed/trim_include/*trim.fq .

nano AlignReads_mcav_include.sh

#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="Align_Mcav_include_out_error"
#SBATCH --output="Align_Mcav_include_out"

module load STAR/2.5.3a-foss-2016b

F=/data/putnamlab/jillashey/Francois_data/Florida/output/STAR/AlignReads_Mcav_include

array1=($(ls $F/*trim.fq))
for i in ${array1[@]}
do
STAR --runMode alignReads --quantMode TranscriptomeSAM --outTmpDir ${i}_TMP --readFilesIn ${i} --genomeDir /data/putnamlab/jillashey/Francois_data/Florida/output/STAR/GenomeIndex_Mcav --twopassMode Basic --twopass1readsN -1 --outStd Log BAM_Unsorted BAM_Quant --outSAMtype BAM Unsorted SortedByCoordinate --outReadsUnmapped Fastx --outFileNamePrefix ${i}.
done 

sbatch AlignReads_mcav_include.sh 
```

Submitted batch job 1717855


### 7) Perform gene counts with stringTie 

####A. cerv

I performed stringTie above with the first batch of samples. I'm going to perform stringTie again with the first and second batch of samples. I moved the results / scripts from the first batch to directories labelled 'old'. 

a) Move BAM files to stringTie folder 

```
cd /data/putnamlab/jillashey/Francois_data/Florida/output/STAR/AlignReads_Acerv_include 
mv *Aligned.sortedByCoord.out.bam ../../../stringTie/Acerv/BAM
```

b) Assemble and estimate reads 

```
cd stringTie/Acerv/BAM

nano stringTie_acerv_assemble.sh

#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="assemble_acerv_out_error"
#SBATCH --output="assemble_acerv_out"

# Running stringTie with all FL samples 

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/Florida/stringTie/Acerv/BAM/

array1=($(ls $F/*bam))
for i in ${array1[@]}; do
	stringtie -G /data/putnamlab/jillashey/genome/Acerv/Acerv_assembly_v1.0.gff3 -e -o ${i}.gtf ${i}
	echo "${i}"
done

sbatch stringTie_acerv_assemble.sh

```
Submitted batch job 1718818

c) Merge stringTie gtf results 

Usually, I do stringTie with all the samples with each annotation file. Now, I think I'm going to move forward only with samples that are assigned to each gff. So only Acerv samples from now on

```
mv 19_T33_Ac_WK.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf 24_T12_Ac_FM.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf 25_ctl1_Ac_GF_1.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf 27_ctl2_Ac_YG_1.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf 31_T22_Ac_UV.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf 35_T43_Ac_MT.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf 37_T13_Ac_ML.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf 38_T23_Ac_IN.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf 41_ctl3_Ac_RN_1.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf 45_T41_Ac_SC_1.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf 47_T31_Ac_JB.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf 52_T11_Ac_II.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf 53_T21_Ac_NH.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf 54_T42_Ac_JQ.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf 57_T32_Ac_NM.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf ../GTF

ls *gtf > acerv_mergelist.txt
cat acerv_mergelist.txt

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

stringtie --merge -p 8 -G /data/putnamlab/jillashey/genome/Acerv/Acerv_assembly_v1.0.gff3 -o stringtie_acerv_merged.gtf acerv_mergelist.txt
```

d) Assess assembly quality

```
module load gffcompare/0.11.5-foss-2018b

gffcompare -r /data/putnamlab/jillashey/genome/Acerv/Acerv_assembly_v1.0.gff3 -o Acerv.merged stringtie_acerv_merged.gtf

48478 reference transcripts loaded.
  48478 query transfrags loaded.
```

e) Re-estimate assembly 

```
nano stringTie_acerv_re-assemble.sh

#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="re-assemble_acerv_out_error"
#SBATCH --output="re-assemble_acerv_out"

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/Florida/stringTie/Acerv/BAM/

array1=($(ls $F/*bam))
for i in ${array1[@]}
do
stringtie -e -G /data/putnamlab/jillashey/genome/Acerv/Acerv_assembly_v1.0.gff3 -o ${i}.merge.gtf ${i}
echo "${i}"
done

sbatch stringTie_acerv_re-assemble.sh

mv 19_T33_Ac_WK.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf 24_T12_Ac_FM.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf 25_ctl1_Ac_GF_1.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf 27_ctl2_Ac_YG_1.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf 31_T22_Ac_UV.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf 35_T43_Ac_MT.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf 37_T13_Ac_ML.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf 38_T23_Ac_IN.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf 41_ctl3_Ac_RN_1.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf 45_T41_Ac_SC_1.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf 47_T31_Ac_JB.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf 52_T11_Ac_II.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf 53_T21_Ac_NH.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf 54_T42_Ac_JQ.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf 57_T32_Ac_NM.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf ../GTF_merge

```
Submitted batch job 1718830

f) Create gene matrix

```
module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/Florida/stringTie/Acerv/GTF_merge/

array2=($(ls *merge.gtf))

for i in ${array2[@]}
do
echo "${i} $F${i}" >> sample_list_acerv.txt
done

python prepDE.py -g gene_count_acerv_only_matrix.csv -i sample_list_acerv.txt
```

g) Secure-copy gene counts onto local computer

```
scp jillashey@bluewaves.uri.edu:/data/putnamlab/jillashey/Francois_data/Florida/stringTie/Acerv/GTF_merge/gene_count_acerv_only_matrix.csv /Users/jillashey/Desktop/
```

#### M. cav

I performed stringTie above with the first batch of samples. I'm going to perform stringTie again with the first and second batch of samples. I moved the results / scripts from the first batch to directories labelled 'old'. 

a) Move BAM files to stringTie folder 

```
cd /data/putnamlab/jillashey/Francois_data/Florida/output/STAR/AlignReads_Mcav_include
mv *Aligned.sortedByCoord.out.bam ../../../stringTie/Mcav/BAM
```

b) Assemble and estimate reads 

```
cd stringTie/Mcav/BAM

nano stringTie_mcav_assemble.sh

#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="assemble_mcav_out_error"
#SBATCH --output="assemble_mcav_out"

# Running stringTie with all FL samples 

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/Florida/stringTie/Mcav/BAM/

array1=($(ls $F/*bam))
for i in ${array1[@]}; do
	stringtie -G /data/putnamlab/jillashey/genome/Mcav/Mcavernosa_annotation/Mcav.gff.annotations.fixed_transcript.gff3 -e -o ${i}.gtf ${i}
	echo "${i}"
done

sbatch stringTie_mcav_assemble.sh

```
Submitted batch job 1718819

c) Merge stringTie gtf results 

Usually, I do stringTie with all the samples with each annotation file. Now, I think I'm going to move forward only with samples that are assigned to each gff. So only Mcav samples from now on

```
mv 20_T12_Mc_PWC.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf 21_T33_Mc_EOU.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf 22_ctl2_Mc_TWF_1.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf 28_ctl1_Mc_GBM_1.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf 29_T23_Mc_PND.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf 33_T43_Mc_RFV.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf 34_T22_Mc_SVS.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf 39_T13_Mc_FJE.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf 42_ctl3_Mc_MGR_1.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf 46_T41_Mc_QYH_1.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf 49_T31_Mc_SWQ.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf 55_T32_Mc_TWP.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf 56_T42_Mc_JAW.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf 58_T21_Mc_EAH.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf 61_T11_Mc_RAP.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf ../GTF


ls *gtf > mcav_mergelist.txt
cat mcav_mergelist.txt

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

stringtie --merge -p 8 -G /data/putnamlab/jillashey/genome/Mcav/Mcavernosa_annotation/Mcav.gff.annotations.fixed_transcript.gff3 -o stringtie_mcav_merged.gtf mcav_mergelist.txt

# Still got these weird errors. Going to keep moving on...
Warning: invalid start coordinate at line:
###			NA	NA				;transcript_id=
Warning: invalid start coordinate at line:
###			NA	NA				;transcript_id=
Warning: invalid start coordinate at line:
###			NA	NA				;transcript_id=
Warning: invalid start coordinate at line:
###			NA	NA				;transcript_id=
```

d) Assess assembly quality

```
module load gffcompare/0.11.5-foss-2018b

gffcompare -r /data/putnamlab/jillashey/genome/Mcav/Mcavernosa_annotation/Mcav.gff.annotations.fixed_transcript.gff3 -o Mcav.merged stringtie_mcav_merged.gtf

25605 reference transcripts loaded.
  25605 query transfrags loaded.
  
# Still got these weird errors. Going to keep moving on...
Warning: invalid start coordinate at line:
###			NA	NA				;transcript_id=
Warning: invalid start coordinate at line:
###			NA	NA				;transcript_id=
Warning: invalid start coordinate at line:
###			NA	NA				;transcript_id=
Warning: invalid start coordinate at line:
###			NA	NA				;transcript_id=
Warning: invalid start coordinate at line:
###			NA	NA				;transcript_id=

```

e) Re-estimate assembly 

```
nano stringTie_mcav_re-assemble.sh

#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="re-assemble_mcav_out_error"
#SBATCH --output="re-assemble_mcav_out"

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/Florida/stringTie/Mcav/BAM/

array1=($(ls $F/*bam))
for i in ${array1[@]}
do
stringtie -e -G /data/putnamlab/jillashey/genome/Mcav/Mcavernosa_annotation/Mcav.gff.annotations.fixed_transcript.gff3 -o ${i}.merge.gtf ${i}
echo "${i}"
done

sbatch stringTie_mcav_re-assemble.sh

mv 20_T12_Mc_PWC.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf 21_T33_Mc_EOU.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf 22_ctl2_Mc_TWF_1.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf 28_ctl1_Mc_GBM_1.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf 29_T23_Mc_PND.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf 33_T43_Mc_RFV.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf 34_T22_Mc_SVS.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf 39_T13_Mc_FJE.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf 42_ctl3_Mc_MGR_1.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf 46_T41_Mc_QYH_1.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf 49_T31_Mc_SWQ.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf 55_T32_Mc_TWP.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf 56_T42_Mc_JAW.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf 58_T21_Mc_EAH.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf 61_T11_Mc_RAP.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf ../GTF_merge
```
Submitted batch job 1718834

f) Create gene matrix

```
module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/Florida/stringTie/Mcav/GTF_merge/

array2=($(ls *merge.gtf))

for i in ${array2[@]}
do
echo "${i} $F${i}" >> sample_list_mcav.txt
done

python prepDE.py -g gene_count_mcav_only_matrix.csv -i sample_list_mcav.txt
```

g) Secure-copy gene counts onto local computer

```
scp jillashey@bluewaves.uri.edu:/data/putnamlab/jillashey/Francois_data/Florida/stringTie/Mcav/GTF_merge/gene_count_mcav_only_matrix.csv /Users/jillashey/Desktop/
```

#### O. fav

I performed stringTie above with the first batch of samples. I'm going to perform stringTie again with the first and second batch of samples. I moved the results / scripts from the first batch to directories labelled 'old'. 

a) Move BAM files to stringTie folder 

```
cd /data/putnamlab/jillashey/Francois_data/Florida/output/STAR/AlignReads_Ofav_include
mv *Aligned.sortedByCoord.out.bam ../../../stringTie/Ofav/BAM
```

b) Assemble and estimate reads 

```
cd stringTie/Ofav/BAM

nano stringTie_ofav_assemble.sh

#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="assemble_ofav_out_error"
#SBATCH --output="assemble_ofav_out"

# Running stringTie with all FL samples 

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/Florida/stringTie/Ofav/BAM/

array1=($(ls $F/*bam))
for i in ${array1[@]}; do
	stringtie -G /data/putnamlab/jillashey/genome/Ofav/GCF_002042975.1_ofav_dov_v1_genomic.gff -e -o ${i}.gtf ${i}
	echo "${i}"
done

sbatch stringTie_ofav_assemble.sh
```

Submitted batch job 1718820


c) Merge stringTie gtf results 

```
mv 17_ctl2_Of_ZTH_1.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf 18_T33_Of_VLL.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf 23_ctl1_Of_CTX_1.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf 26_T12_Of_WCL.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf 30_T23_Of_RPG.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf 32_T22_Of_EVR.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf 36_T43_Of_JJN.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf 40_T13_Of_GWS.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf 43_ctl3_Of_JVP_1.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf 44_T41_Of_PVT_1.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf 48_T31_Of_JNO.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf 50_T21_Of_YZB.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf 51_T42_Of_UOF.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf 59_T11_Of_TQP.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf 60_T32_Of_WXY.fastq.trim.fq.Aligned.sortedByCoord.out.bam.gtf ../GTF

ls *gtf > ofav_mergelist.txt
cat ofav_mergelist.txt

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

stringtie --merge -p 8 -G /data/putnamlab/jillashey/genome/Ofav/GCF_002042975.1_ofav_dov_v1_genomic.gff -o stringtie_ofav_merged.gtf ofav_mergelist.txt
```

d) Assess assembly quality

```
module load gffcompare/0.11.5-foss-2018b

gffcompare -r /data/putnamlab/jillashey/genome/Ofav/GCF_002042975.1_ofav_dov_v1_genomic.gff -o Ofav.merged stringtie_ofav_merged.gtf

 37786 reference transcripts loaded.
  5 duplicate reference transcripts discarded.
  37782 query transfrags loaded.
```

e) Re-estimate assembly 

```
nano stringTie_ofav_re-assemble.sh

#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="re-assemble_ofav_out_error"
#SBATCH --output="re-assemble_ofav_out"

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/Florida/stringTie/Ofav/BAM/

array1=($(ls $F/*bam))
for i in ${array1[@]}
do
stringtie -e -G /data/putnamlab/jillashey/genome/Ofav/GCF_002042975.1_ofav_dov_v1_genomic.gff -o ${i}.merge.gtf ${i}
echo "${i}"
done

sbatch stringTie_ofav_re-assemble.sh

mv 17_ctl2_Of_ZTH_1.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf 18_T33_Of_VLL.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf 23_ctl1_Of_CTX_1.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf 26_T12_Of_WCL.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf 30_T23_Of_RPG.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf 32_T22_Of_EVR.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf 36_T43_Of_JJN.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf 40_T13_Of_GWS.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf 43_ctl3_Of_JVP_1.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf 44_T41_Of_PVT_1.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf 48_T31_Of_JNO.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf 50_T21_Of_YZB.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf 51_T42_Of_UOF.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf 59_T11_Of_TQP.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf 60_T32_Of_WXY.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf ../GTF_merge
```
Submitted batch job 1718849

f) Create gene matrix

```
module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/Florida/stringTie/Ofav/GTF_merge/

array2=($(ls *merge.gtf))

for i in ${array2[@]}
do
echo "${i} $F${i}" >> sample_list_ofav.txt
done

python prepDE.py -g gene_count_ofav_only_matrix.csv -i sample_list_ofav.txt
```

g) Secure-copy gene counts onto local computer

```
scp jillashey@bluewaves.uri.edu:/data/putnamlab/jillashey/Francois_data/Florida/stringTie/Ofav/GTF_merge/gene_count_ofav_only_matrix.csv /Users/jillashey/Desktop/
```

