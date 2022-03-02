20220228

Going through all data files and code on bluewaves for FL species 

```
cd /data/putnamlab/jillashey/Francois_data/Florida
ls
blast  data  fastqc_results  multiqc_results  output  scripts  stringTie
```

Species: _Acropora cervicornis_, _Montestraea cavernosa_, and _Orbicella faveolata_

```
cd blast 
ls
acerv  extract_proteins.py  mcav  mcav_transcript_remote_blast_maxseqs_out  mcav_transcript_remote_blast_maxseqs_out_error  pdam  plob
```

Removed entire blast folder, old results

```
cd data 
ls
Illumina_adapter_reads_PE_SE.fa  raw  sampleNames  slurm-1690400.out  slurm-1690401.out  trimmed
```

Keep adapter reads, remove s*

Check transfer 

```
nano check_transfer.sh

#!/bin/bash
#SBATCH -t 18:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@my.uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="check_transfer_out_error"
#SBATCH --output="check_transfer_out"

md5sum /data/putnamlab/jillashey/Francois_data/Florida/data/raw/*fastq > checkmd5.md5

sbatch check_transfer.sh 
```








QC

Raw 

```
module load MultiQC/1.7-foss-2018b-Python-2.7.15
[jillashey@bluewaves multiqc_results]$ multiqc /data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/raw/*fastqc.zip -o /data/putnamlab/jillashey/Francois_data/Florida/multiqc_results/
[WARNING]         multiqc : MultiQC Version v1.12 now available!
[INFO   ]         multiqc : This is MultiQC v1.7
[INFO   ]         multiqc : Template    : default
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/raw/17_ctl2_Of_ZTH_1_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/raw/18_T33_Of_VLL_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/raw/19_T33_Ac_WK_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/raw/20_T12_Mc_PWC_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/raw/21_T33_Mc_EOU_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/raw/22_ctl2_Mc_TWF_1_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/raw/23_ctl1_Of_CTX_1_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/raw/24_T12_Ac_FM_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/raw/25_ctl1_Ac_GF_1_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/raw/26_T12_Of_WCL_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/raw/27_ctl2_Ac_YG_1_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/raw/28_ctl1_Mc_GBM_1_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/raw/29_T23_Mc_PND_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/raw/30_T23_Of_RPG_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/raw/31_T22_Ac_UV_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/raw/32_T22_Of_EVR_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/raw/33_T43_Mc_RFV_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/raw/34_T22_Mc_SVS_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/raw/35_T43_Ac_MT_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/raw/36_T43_Of_JJN_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/raw/37_T13_Ac_ML_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/raw/38_T23_Ac_IN_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/raw/39_T13_Mc_FJE_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/raw/40_T13_Of_GWS_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/raw/41_ctl3_Ac_RN_1_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/raw/42_ctl3_Mc_MGR_1_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/raw/43_ctl3_Of_JVP_1_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/raw/44_T41_Of_PVT_1_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/raw/45_T41_Ac_SC_1_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/raw/46_T41_Mc_QYH_1_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/raw/47_T31_Ac_JB_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/raw/48_T31_Of_JNO_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/raw/49_T31_Mc_SWQ_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/raw/50_T21_Of_YZB_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/raw/51_T42_Of_UOF_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/raw/52_T11_Ac_II_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/raw/53_T21_Ac_NH_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/raw/54_T42_Ac_JQ_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/raw/55_T32_Mc_TWP_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/raw/56_T42_Mc_JAW_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/raw/57_T32_Ac_NM_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/raw/58_T21_Mc_EAH_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/raw/59_T11_Of_TQP_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/raw/60_T32_Of_WXY_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/raw/61_T11_Mc_RAP_fastqc.zip'
[INFO   ]          fastqc : Found 45 reports
[INFO   ]         multiqc : Compressing plot data
[INFO   ]         multiqc : Report      : multiqc_report.html
[INFO   ]         multiqc : Data        : multiqc_data
[INFO   ]         multiqc : MultiQC complete
```

trimmed 

```
multiqc /data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/trimmed/*fastqc.zip -o /data/putnamlab/jillashey/Francois_data/Florida/multiqc_results/
[WARNING]         multiqc : MultiQC Version v1.12 now available!
[INFO   ]         multiqc : This is MultiQC v1.7
[INFO   ]         multiqc : Template    : default
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/trimmed/17_ctl2_Of_ZTH_1.fastq.trim_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/trimmed/18_T33_Of_VLL.fastq.trim_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/trimmed/19_T33_Ac_WK.fastq.trim_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/trimmed/20_T12_Mc_PWC.fastq.trim_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/trimmed/21_T33_Mc_EOU.fastq.trim_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/trimmed/22_ctl2_Mc_TWF_1.fastq.trim_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/trimmed/23_ctl1_Of_CTX_1.fastq.trim_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/trimmed/24_T12_Ac_FM.fastq.trim_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/trimmed/25_ctl1_Ac_GF_1.fastq.trim_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/trimmed/26_T12_Of_WCL.fastq.trim_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/trimmed/27_ctl2_Ac_YG_1.fastq.trim_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/trimmed/28_ctl1_Mc_GBM_1.fastq.trim_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/trimmed/29_T23_Mc_PND.fastq.trim_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/trimmed/30_T23_Of_RPG.fastq.trim_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/trimmed/31_T22_Ac_UV.fastq.trim_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/trimmed/32_T22_Of_EVR.fastq.trim_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/trimmed/33_T43_Mc_RFV.fastq.trim_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/trimmed/34_T22_Mc_SVS.fastq.trim_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/trimmed/35_T43_Ac_MT.fastq.trim_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/trimmed/36_T43_Of_JJN.fastq.trim_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/trimmed/37_T13_Ac_ML.fastq.trim_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/trimmed/38_T23_Ac_IN.fastq.trim_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/trimmed/39_T13_Mc_FJE.fastq.trim_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/trimmed/40_T13_Of_GWS.fastq.trim_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/trimmed/41_ctl3_Ac_RN_1.fastq.trim_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/trimmed/42_ctl3_Mc_MGR_1.fastq.trim_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/trimmed/43_ctl3_Of_JVP_1.fastq.trim_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/trimmed/44_T41_Of_PVT_1.fastq.trim_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/trimmed/45_T41_Ac_SC_1.fastq.trim_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/trimmed/46_T41_Mc_QYH_1.fastq.trim_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/trimmed/47_T31_Ac_JB.fastq.trim_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/trimmed/48_T31_Of_JNO.fastq.trim_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/trimmed/49_T31_Mc_SWQ.fastq.trim_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/trimmed/50_T21_Of_YZB.fastq.trim_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/trimmed/51_T42_Of_UOF.fastq.trim_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/trimmed/52_T11_Ac_II.fastq.trim_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/trimmed/53_T21_Ac_NH.fastq.trim_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/trimmed/54_T42_Ac_JQ.fastq.trim_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/trimmed/55_T32_Mc_TWP.fastq.trim_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/trimmed/56_T42_Mc_JAW.fastq.trim_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/trimmed/57_T32_Ac_NM.fastq.trim_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/trimmed/58_T21_Mc_EAH.fastq.trim_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/trimmed/59_T11_Of_TQP.fastq.trim_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/trimmed/60_T32_Of_WXY.fastq.trim_fastqc.zip'
[INFO   ]         multiqc : Searching '/data/putnamlab/jillashey/Francois_data/Florida/fastqc_results/trimmed/61_T11_Mc_RAP.fastq.trim_fastqc.zip'
[INFO   ]          fastqc : Found 45 reports
[INFO   ]         multiqc : Compressing plot data
[INFO   ]         multiqc : Report      : multiqc_report.html
[INFO   ]         multiqc : Data        : multiqc_data
[INFO   ]         multiqc : MultiQC complete
```