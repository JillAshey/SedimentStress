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

![](https://raw.githubusercontent.com/JillAshey/SedimentStress/master/Images/QC/trimmed/fastqc_sequence_counts_plot.png?token=APHKO32IMCD4U5AQOXNGYMS7J726W)

![](https://raw.githubusercontent.com/JillAshey/SedimentStress/master/Images/QC/trimmed/fastqc_per_sequence_gc_content_plot.png?token=APHKO37KY2NCDO56RVLUQPC7J73CC)

![](https://raw.githubusercontent.com/JillAshey/SedimentStress/master/Images/QC/trimmed/fastqc_per_base_sequence_quality_plot.png?token=APHKO32JTY6U2JJXCY2J3U27J73JW)

![](https://raw.githubusercontent.com/JillAshey/SedimentStress/master/Images/QC/trimmed/fastqc_overrepresented_sequencesi_plot.png?token=APHKO33VGTBNNZAGIKCKSD27J73OY)

![](https://raw.githubusercontent.com/JillAshey/SedimentStress/master/Images/QC/trimmed/fastqc_adapter_content_plot.png?token=APHKO334NFOXWGH7KFVQZRK7J73PY)

Adapter content looks kinda weird...maybe because Francois may have already trimmed some of the samples? But not sure about that yet

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
Submitted batch job 1691569
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

Might have to do some gff fixing for this gff file

```
module load STAR/2.5.3a-foss-2016b

STAR --runThreadN 10 --runMode genomeGenerate --genomeDir /data/putnamlab/jillashey/Francois_data/Florida/output/STAR/GenomeIndex_Mcav --genomeFastaFiles /data/putnamlab/jillashey/genome/Mcav/Mcavernosa_July2018.fasta --sjdbGTFfile /data/putnamlab/jillashey/genome/Mcav/Mcavernosa_annotation/Mcavernosa.maker.coding.gff3

```
Error: terminate called after throwing an instance of 'std::out_of_range'
  what():  vector::_M_range_check: __n (which is 0) >= this->size() (which is 0)
Aborted

Need to go fix gff in r. 

Added transcript_id to gene col. Should work now

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

I am interested to trying --quantMode GeneCounts in the STAR align reads step. Going to try with an Ofav sample and see what the output is.

```
mkdir test_quantMode
cp 32_T22_Of_EVR.fastq.trim.fq test_quantMode
cd test_quantMode

nano AlignReads_ofav_32_test.sh

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

STAR --runMode alignReads --quantMode GeneCounts TranscriptomeSAM --outTmpDir 32_T22_Of_EVR.fastq.trim.fq_TMP --readFilesIn 32_T22_Of_EVR.fastq.trim.fq --genomeDir /data/putnamlab/jillashey/Francois_data/Florida/output/STAR/GenomeIndex_Ofav --twopassMode Basic --twopass1readsN -1 --outStd Log BAM_Unsorted BAM_Quant --outSAMtype BAM Unsorted SortedByCoordinate --outReadsUnmapped Fastx --outFileNamePrefix 32_T22_Of_EVR.fastq.trim.fq.

sbatch AlignReads_ofav_32_test.sh 
```
Submitted batch job 1692813

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


### 6-repeat) Align reads with STAR - Mcav only

a) Generate genome index

```
module load STAR/2.5.3a-foss-2016b

STAR --runThreadN 10 --runMode genomeGenerate --genomeDir /data/putnamlab/jillashey/Francois_data/Florida/output/STAR/GenomeIndex_Mcav_gene.diff.Pos --genomeFastaFiles /data/putnamlab/jillashey/genome/Mcav/Mcavernosa_July2018.fasta --sjdbGTFfile /data/putnamlab/jillashey/genome/Mcav/Mcavernosa_annotation/Mcav.gff.annotations.fixed_gene.diff.Pos.gff3

```
Fatal INPUT FILE error, no exon lines in the GTF file: /data/putnamlab/jillashey/genome/Mcav/Mcavernosa_annotation/Mcav.gff.annotations.fixed_gene.diff.Pos.gff3
Solution: check the formatting of the GTF file, it must contain some lines with exon in the 3rd column.
          Make sure the GTF file is unzipped.
          If exons are marked with a different word, use --sjdbGTFfeatureExon .

But the gff has exons...........confused. Going to use the --sjdbGTFfeatureExon .


```
module load STAR/2.5.3a-foss-2016b

STAR --runThreadN 10 --runMode genomeGenerate --genomeDir /data/putnamlab/jillashey/Francois_data/Florida/output/STAR/GenomeIndex_Mcav_gene.diff.Pos --genomeFastaFiles /data/putnamlab/jillashey/genome/Mcav/Mcavernosa_July2018.fasta --sjdbGTFfile /data/putnamlab/jillashey/genome/Mcav/Mcavernosa_annotation/Mcav.gff.annotations.fixed_gene.diff.Pos.gff3 --sjdbGTFfeatureExon .

```

same error.......idk












b) Align reads to genome

```
mkdir AlignReads_Mcav_gene.diff.Pos
cd AlignReads_Mcav_gene.diff.Pos
ln -s /data/putnamlab/jillashey/Francois_data/Florida/data/trimmed/*trim.fq .

nano AlignReads_gene.diff.Pos_mcav.sh

#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="Align_gene.diff.Pos_mav_out_error"
#SBATCH --output="Align_gene.diff.Pos_mcav_out"

module load STAR/2.5.3a-foss-2016b

F=/data/putnamlab/jillashey/Francois_data/Florida/output/STAR/AlignReads_Mcav_gene.diff.Pos

array1=($(ls $F/*trim.fq))
for i in ${array1[@]}
do
STAR --runMode alignReads --quantMode TranscriptomeSAM --outTmpDir ${i}_TMP --readFilesIn ${i} --genomeDir /data/putnamlab/jillashey/Francois_data/Florida/output/STAR/GenomeIndex_Mcav_gene.diff.Pos --twopassMode Basic --twopass1readsN -1 --outStd Log BAM_Unsorted BAM_Quant --outSAMtype BAM Unsorted SortedByCoordinate --outReadsUnmapped Fastx --outFileNamePrefix ${i}.
done 

sbatch AlignReads_gene.diff.Pos_mcav.sh 
```

Fatal INPUT FILE error, no exon lines in the GTF file: /data/putnamlab/jillashey/genome/Mcav/Mcavernosa_annotation/Mcav.gff.annotations.fixed_gene.diff.Pos.gff3
Solution: check the formatting of the GTF file, it must contain some lines with exon in the 3rd column.
          Make sure the GTF file is unzipped.
          If exons are marked with a different word, use --sjdbGTFfeatureExon .

But the gff has exons...........confused. Going to use the --sjdbGTFfeatureExon .
















For now, continuing on with gff that still has - gene lengths 

d) Assess assembly quality

```
module load gffcompare/0.11.5-foss-2018b

gffcompare -r /data/putnamlab/jillashey/genome/Mcav/Mcavernosa_annotation/Mcav.gff.annotations.fixed_transcript.gff3 -o Mcav.merged stringtie_mcav_merged.gtf


Warning: invalid start coordinate at line:
###			NA	NA				;transcript_id=
Warning: invalid start coordinate at line:
###			NA	NA				;transcript_id=
 25605 reference transcripts loaded.
  25605 query transfrags loaded.
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

mv *merge.gtf ../GTF_merge
```

Submitted batch job 1693072


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

python prepDE.py -g gene_count_mcav_matrix.csv -i sample_list_mcav.txt
```

g) Secure-copy gene counts onto local computer

```
scp jillashey@bluewaves.uri.edu:/data/putnamlab/jillashey/Francois_data/Florida/stringTie/Mcav/GTF_merge/gene_count_mcav_matrix.csv /Users/jillashey/Desktop/Putnamlab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/star/
```

