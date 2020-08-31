## Pipeline for FL species - sediment stress 

#### 1) Check file integrity 

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

#### 2) Run FastQC

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

#### 3) Run MultiQC

a) Make folders for raw MultiWC results

```
cd Francois_data/Florida
mkdir multiqc_results/raw
```

b) Run MultiQC. Pretty fast, so don't need to submit job for it 

```
module load MultiQC/1.7-foss-2018b-Python-2.7.15
multiqc /data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/raw/fastqc.zip -o /data/putnamlab/jillashey/Francois_data/Florida/multiqc_results/raw
```

c) Copy MultiQC files to local computer

```
scp -r jillashey@bluewaves.uri.edu:/data/putnamlab/jillashey/Francois_data/Florida/multiqc_results/raw/* 

```
