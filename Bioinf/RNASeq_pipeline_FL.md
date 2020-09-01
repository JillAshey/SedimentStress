## Pipeline for FL species - sediment stress 

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
18_T33_Of_VLL.txt.gz:7191417
19_T33_Ac_WK.txt.gz:19665420
20_T12_Mc_PWC.txt.gz:18018163
21_T33_Mc_EOU.txt.gz:15475674
24_T12_Ac_FM.txt.gz:13786799
26_T12_Of_WCL.txt.gz:13436090
29_T23_Mc_PND.txt.gz:0
30_T23_Of_RPG.txt.gz:0
31_T22_Ac_UV.txt.gz:0
32_T22_Of_EVR.txt.gz:0
34_T22_Mc_SVS.txt.gz:12609087
35_T43_Ac_MT.txt.gz:23846590
36_T43_Of_JJN.txt.gz:13551094
37_T13_Ac_ML.txt.gz:14171723
38_T23_Ac_IN.txt.gz:21348241
39_T13_Mc_FJE.txt.gz:16994083
40_T13_Of_GWS.txt.gz:14510120
gzip: 43_ctl3_Of_JVP_2.txt.gz: unexpected end of file
43_ctl3_Of_JVP_2.txt.gz:11461086
47_T31_Ac_JB.txt.gz:14229140
48_T31_Of_JNO.txt.gz:12280565
49_T31_Mc_SWQ.txt.gz:14934691
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
29_T23_Mc_PND.txt.gz:15499497
30_T23_Of_RPG.txt.gz:8218881
31_T22_Ac_UV.txt.gz:13317263
32_T22_Of_EVR.txt.gz:12667861
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
50_T21_Of_YZB.txt.gz:15084847
51_T42_Of_UOF.txt.gz:14442222
52_T11_Ac_II.txt.gz:17303453
53_T21_Ac_NH.txt.gz:10306635
54_T42_Ac_JQ.txt.gz:23542186
55_T32_Mc_TWP.txt.gz:19313398
56_T42_Mc_JAW.txt.gz:11743795
57_T32_Ac_NM.txt.gz:11679142
58_T21_Mc_EAH.txt.gz:20792838
59_T11_Of_TQP.txt.gz:19297200
60_T32_Of_WXY.txt.gz:19349770
61_T11_Mc_RAP.txt.gz:16873902
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
```
```
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
```
```
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

![](https://raw.githubusercontent.com/JillAshey/SedimentStress/master/Images/QC/raw/FL/fastqc_sequence_counts_plot.png?token=APHKO32WXRFGOWEOOB7GK4S7JXJTU)

![](https://raw.githubusercontent.com/JillAshey/SedimentStress/master/Images/QC/raw/FL/fastqc_per_sequence_gc_content_plot.png?token=APHKO3ZF2KKDAHO3YOSUG4S7JXJRG)

![](https://raw.githubusercontent.com/JillAshey/SedimentStress/master/Images/QC/raw/FL/fastqc_per_base_sequence_quality_plot.png?token=APHKO33J7GUAQSUGQ75LGW27JXJO6)

![](https://raw.githubusercontent.com/JillAshey/SedimentStress/master/Images/QC/raw/FL/fastqc_overrepresented_sequencesi_plot.png?token=APHKO33KBAXH42LTQTUCNCC7JXJNS)

![](https://raw.githubusercontent.com/JillAshey/SedimentStress/master/Images/QC/raw/FL/fastqc_adapter_content_plot.png?token=APHKO32ATSLZDK5WULLAVYK7JXJHM)

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
```

```
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
```
```
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
```
Submitted batch job 1690400 for HISEQ and 1690401 for HWI

c) Run FastQC on trimmed data

```
nano fastqc_trimmed.sh
```
```
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
```
```
sbatch fastqc_trimmed.sh
```
Submitted batch job 1690402

d) Run MultiQC on trimmed data and moved files 

```
module load MultiQC/1.7-foss-2018b-Python-2.7.15
multiqc /data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/trimmed/*fastqc.zip -o /data/putnamlab/jillashey/Francois_data/Florida/multiqc_results/trimmed

scp -r jillashey@bluewaves.uri.edu:/data/putnamlab/jillashey/Francois_data/Florida/multiqc_results/trimmed/multiqc_data /Users/jillashey/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/QC/trimmed

scp jillashey@bluewaves.uri.edu:/data/putnamlab/jillashey/Francois_data/Florida/multiqc_results/trimmed/multiqc_report.html /Users/jillashey/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/QC/trimmed
```

### 6) Align reads with STAR

#### O. fav

a) Generate genome index

```
module load STAR/2.5.3a-foss-2016b

STAR --runThreadN 10 --runMode genomeGenerate --genomeDir /data/putnamlab/jillashey/Francois_data/Florida/output/STAR/GenomeIndex_Ofav --genomeFastaFiles /data/putnamlab/jillashey/genome/Ofav/GCF_002042975.1_ofav_dov_v1_genomic.fna --sjdbGTFfile /data/putnamlab/jillashey/genome/Ofav/GCF_002042975.1_ofav_dov_v1_genomic.gff
```

b) Align reads to genome

Test with one sample first 

```
nano staralign_18_T33_Of_VLL_test.sh

#!/bin/bash
#SBATCH -t 48:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="AlignReads_18_T33_Of_VLL_test_out_error"
#SBATCH --output="AlignReads_18_T33_Of_VLL_test_out"

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads --quantMode TranscriptomeSAM --outTmpDir /data/putnamlab/jillashey/Francois_data/Florida/output/STAR/AlignReads_Ofav/18_T33_Of_VLL.fastq.trim.fq_TMP --readFilesIn /data/putnamlab/jillashey/Francois_data/Florida/data/trimmed/18_T33_Of_VLL.fastq.trim.fq --genomeDir /data/putnamlab/jillashey/Francois_data/Florida/output/STAR/GenomeIndex_Ofav --twopassMode Basic --twopass1readsN -1 --outStd Log BAM_Unsorted BAM_Quant --outSAMtype BAM Unsorted SortedByCoordinate --outReadsUnmapped Fastx --outFileNamePrefix /data/putnamlab/jillashey/Francois_data/Florida/output/STAR/AlignReads_Ofav/18_T33_Of_VLL.fastq.trim.fq

sbatch staralign_18_T33_Of_VLL_test.sh
```

```
nano staralign_Ofav.sh

#!/bin/bash
#SBATCH -t 48:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="GenomeIndex_Ofav_out_error"
#SBATCH --output="GenomeIndex_Ofac_out"

module load STAR/2.5.3a-foss-2016b

FILENAME=/data/putnamlab/jillashey/Francois_data/Florida/data/trimmed/sampleNames

STAR --runMode alignReads --quantMode TranscriptomeSAM --outTmpDir /data/putnamlab/jillashey/Francois_data/Florida/output/STAR/${FILENAME}_TMP --readFilesIn /data/putnamlab/jillashey/Francois_data/Florida/data/trimmed/${FILENAME} --genomeDir /data/putnamlab/jillashey/Francois_data/Florida/output/STAR/GenomeIndex_Ofav --twopassMode Basic --twopass1readsN -1 --outStd Log BAM_Unsorted BAM_Quant --outSAMtype BAM Unsorted SortedByCoordinate --outReadsUnmapped Fastx --outFileNamePrefix /data/putnamlab/jillashey/Francois_data/Florida/output/STAR/AlignReads_Ofav/${FILENAME}.

sbatch staralign_Ofav.sh
```

Submitted batch job 1690404 not working




to isolate name

```
ls -l *trim.fq | awk -F. '{print $9}' > test

awk '{print $9}' test > sampleNames

ls -l /data/putnamlab/jillashey/Francois_data/Florida/data/trimmed/*trim.fq | awk '{print $9}' | awk -F/ '{print $1}' > blah

ls -l /data/putnamlab/jillashey/Francois_data/Florida/data/trimmed/*trim.fq | awk -F/ '{print $9}' > blah

```

Doing filenames with array

```
# To isolate sample names: 
ls -l /data/putnamlab/jillashey/Francois_data/Florida/data/trimmed/*trim.fq | awk -F/ '{print $9}' > sampleNames


nano staralign_Ofav.sh

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

array1=sampleNames

for i in ${array1[@]}
do
STAR --runMode alignReads --quantMode TranscriptomeSAM --outTmpDir /data/putnamlab/jillashey/Francois_data/Florida/output/STAR/${i}_TMP --readFilesIn /data/putnamlab/jillashey/Francois_data/Florida/data/trimmed/${i} --genomeDir /data/putnamlab/jillashey/Francois_data/Florida/output/STAR/GenomeIndex_Ofav --twopassMode Basic --twopass1readsN -1 --outStd Log BAM_Unsorted BAM_Quant --outSAMtype BAM Unsorted SortedByCoordinate --outReadsUnmapped Fastx --outFileNamePrefix /data/putnamlab/jillashey/Francois_data/Florida/output/STAR/AlignReads_Ofav/${i}.
done

sbatch --array 1-3 staralign_Ofav.sh
```
Giving me this error:
ls: cannot access data/putnamlab/jillashey/Francois_data/Florida/data/trimmed/*trim.fq: No such file or directory

```
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

for file in /data/putnamlab/jillashey/Francois_data/Florida/data/trimmed/*trim.fq
do
STAR --runMode alignReads --quantMode TranscriptomeSAM --outTmpDir /data/putnamlab/jillashey/Francois_data/Florida/output/STAR/${file}_TMP --readFilesIn /data/putnamlab/jillashey/Francois_data/Florida/data/trimmed/${file} --genomeDir /data/putnamlab/jillashey/Francois_data/Florida/output/STAR/GenomeIndex_Ofav --twopassMode Basic --twopass1readsN -1 --outStd Log BAM_Unsorted BAM_Quant --outSAMtype BAM Unsorted SortedByCoordinate --outReadsUnmapped Fastx --outFileNamePrefix /data/putnamlab/jillashey/Francois_data/Florida/output/STAR/AlignReads_Ofav/${file}.
done
```
Giving me error that it cannot find files names because it is combining them together in some weird way. Think its something to do with filename input and paths

```
nano Ofav_test.sh

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

array1=blah

for i in ${array1[@]}
do
STAR --runMode alignReads --quantMode TranscriptomeSAM --outTmpDir /data/putnamlab/jillashey/Francois_data/Florida/output/STAR/${i}_TMP --readFilesIn /data/putnamlab/jillashey/Francois_data/Florida/data/trimmed/${i} --genomeDir /data/putnamlab/jillashey/Francois_data/Florida/output/STAR/GenomeIndex_Ofav --twopassMode Basic --twopass1readsN -1 --outStd Log BAM_Unsorted BAM_Quant --outSAMtype BAM Unsorted SortedByCoordinate --outReadsUnmapped Fastx --outFileNamePrefix /data/putnamlab/jillashey/Francois_data/Florida/output/STAR/AlignReads_Ofav/${i}.
done 
```

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
This seems to have worked. Issues with the paths I think

#### A. cerv

a) Generate genome index

```
module load STAR/2.5.3a-foss-2016b

STAR --runThreadN 10 --runMode genomeGenerate --genomeDir /data/putnamlab/jillashey/Francois_data/Florida/output/STAR/GenomeIndex_Acerv --genomeFastaFiles /data/putnamlab/jillashey/genome/Acerv/Acerv_assembly_v1.0_171209.fasta --sjdbGTFfile /data/putnamlab/jillashey/genome/Acerv/Acerv_assembly_v1.0.gff3 

```
Error: terminate called after throwing an instance of 'std::out_of_range'
  what():  vector::_M_range_check: __n (which is 0) >= this->size() (which is 0)
Aborted
Based on previous runs, this is an issue with the ID= in the gff. Gotta fix in R

Added transcript_id to gene col. Should work now

```
module load STAR/2.5.3a-foss-2016b

STAR --runThreadN 10 --runMode genomeGenerate --genomeDir /data/putnamlab/jillashey/Francois_data/Florida/output/STAR/GenomeIndex_Acerv --genomeFastaFiles /data/putnamlab/jillashey/genome/Acerv/Acerv_assembly_v1.0_171209.fasta --sjdbGTFfile /data/putnamlab/jillashey/genome/Acerv/Acerv.GFFannotations.fixed_transcript.gff3

```






b) Align reads to genome

```
mkdir AlignReads_Acerv
cd AlignReads_Acerv
ln -s /data/putnamlab/jillashey/Francois_data/Florida/data/trimmed/*trim.fq .

nano Align_script.sh

#!/bin/bash
#SBATCH -t 100:00:00
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

