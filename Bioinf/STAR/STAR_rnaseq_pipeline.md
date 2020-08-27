# RNA-Seq STAR pipeline

#### Author: Jill Ashey

### 1) Obtain raw data and give it quick check

a) Make directories for data

cd /data/putnamlab/jillashey/
mkdir Francois_data
cd Francois_data
mkdir data/raw 
mkdir data/trimmed

b) Copy Francois's data from Tufts to URI server

```
scp -r jashey@linux.eecs.tufts.edu:/r/corals/Francois/Public/KML_sedim_RNASeq_fastq/Fastq/11698R/*.fastq.gz jillashey@bluewaves.uri.edu:/data/putnamlab/jillashey/Francois_data/data/raw/

# Illumina adapter
scp jashey@linux.eecs.tufts.edu:/r/corals/Francois_data/fasta/Illumina_adapter_reads_PE_SE.fa jillashey@bluewaves.uri.edu:/data/putnamlab/jillashey/Francois_data/data
```

c) Check file count, raw read count, and data integrity 

```
ls -1 | wc -l
```
64 total files

Check number of raw reads per file

```
nano read-count.sh
```

```
#!/bin/bash
#SBATCH -t 4:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@my.uri.edu
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Francois_data/scripts
#SBATCH --error="read-count_raw_error"
#SBATCH --output="read-count_raw"

for file in /data/putnamlab/jillashey/Francois_data/data/raw/*.fastq.gz
do
echo $file >> read-counts.txt
zcat $file | echo $((`wc -l`/4)) >> read-counts.txt
done
```

Submitted batch job 1652373

Looks normal except for /data/putnamlab/jillashey/Francois_data/data/raw/31_2.fastq.gz
6132

```
md5sum *fastq.gz > raw_checksum.md5
md5sum -c raw_checksum.md5
```
Print out 'OK' next to each file

### 2) Check quality with fastqc

Make directories for fastqc results 

```
cd Francois_data
mkdir fastqc_results/raw
```

Write fastqc script

```
cd Francois_data/scripts
nano fastqc_raw.sh
```

```
#!/bin/bash
#SBATCH -t 18:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@my.uri.edu
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Francois_data/scripts
#SBATCH --error="fastqc_out_raw_error"
#SBATCH --output="fastqc_out_raw"

module load FastQC/0.11.8-Java-1.8 #another version in Bluewaves if this one doesn't work

for file in /data/putnamlab/jillashey/Francois_data/data/raw/*.fastq.gz
do
fastqc $file --outdir /data/putnamlab/jillashey/Francois_data/fastqc_results/raw/
done
```

Submitted batch job 1652378

After job is complete, make sure all data is there 

Check number of files 
```
ls -1 | wc -l 
```

2x the number of files, so 128 total files here 

### 3) Run MultiQC

Make directory for multiqc results

```
mkdir Connellydata/data/multiqc_results/raw
```
Load module and run multiqc

```
module load MultiQC/1.7-foss-2018b-Python-2.7.15
multiqc /data/putnamlab/jillashey/Francois_data/fastqc_results/raw/*fastqc.zip
```

Secure-copy multiqc results onto local computer 

```
scp jillashey@bluewaves.uri.edu:/data/putnamlab/jillashey/Francois_data/multiqc/raw/multiqc_report.html /Users/jillashey/Desktop/
scp -r jillashey@bluewaves.uri.edu:/data/putnamlab/jillashey/Francois_data/multiqc/raw/multiqc_data /Users/jillashey/Desktop/MultiQC 
```

![fastqc_sequence_counts_plot.png](https://raw.githubusercontent.com/hputnam/Tufts_URI_CSM_RNASeq/master/_images/SenecaQC/raw/fastqc_sequence_counts_plot.png?token=APHKO34IJAACHJ2AGYU24LK7EIFF2)

![fastqc_per_base_sequence_quality_plot.png](https://raw.githubusercontent.com/hputnam/Tufts_URI_CSM_RNASeq/master/_images/SenecaQC/raw/fastqc_per_base_sequence_quality_plot.png?token=APHKO37RPDSFLLTDRZUU5LK7EIFJQ)

![fastqc_per_sequence_quality_scores_plot.png](https://raw.githubusercontent.com/hputnam/Tufts_URI_CSM_RNASeq/master/_images/SenecaQC/raw/fastqc_per_sequence_quality_scores_plot.png?token=APHKO35CEIDUQRA4CE5WBYC7EIFPG)

![fastqc_per_sequence_base_content_plot.png](https://raw.githubusercontent.com/hputnam/Tufts_URI_CSM_RNASeq/master/_images/SenecaQC/raw/fastqc_per_sequence_base_content_plot.png?token=APHKO32HK3AWC2SYYIUY4PS7EIFR6)

![fastqc_per_sequence_gc_content_plot.png](https://raw.githubusercontent.com/hputnam/Tufts_URI_CSM_RNASeq/master/_images/SenecaQC/raw/fastqc_per_sequence_gc_content_plot.png?token=APHKO34FVCEMH5ELF4F3UNK7EIFZC)

![fastqc_sequence_duplication_levels_plot.png](https://raw.githubusercontent.com/hputnam/Tufts_URI_CSM_RNASeq/master/_images/SenecaQC/raw/fastqc_sequence_duplication_levels_plot.png?token=APHKO36GTRF3MKZNRB4ZFZK7EIF4W)

![fastqc_overrepresented_sequencesi_plot.png](https://raw.githubusercontent.com/hputnam/Tufts_URI_CSM_RNASeq/master/_images/SenecaQC/raw/fastqc_overrepresented_sequencesi_plot.png?token=APHKO36NX7UULKC3PM4HJBK7EIGAM)

![fastqc_adapter_content_plot.png](https://raw.githubusercontent.com/hputnam/Tufts_URI_CSM_RNASeq/master/_images/SenecaQC/raw/fastqc_adapter_content_plot.png?token=APHKO365LPURJG5HKQNUBD27EIGC4)

### 4) Trim reads 

```
#!/bin/bash
#SBATCH -t 48:00:00
#SBATCH --nodes=1 --ntasks-per-node=2
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@my.uri.edu
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Francois_data/scripts
#SBATCH --error="trimmomatic_out_error"
#SBATCH --output="trimmomatic_out"

module load Trimmomatic/0.38-Java-1.8

for file in /data/putnamlab/jillashey/Francois_data/data/raw/*.fastq.gz
do
  	java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.38.jar SE -phred33 $file $file.trim.fq.gz ILLUMINACLIP:/data/putnamlab/jillashey/Francois_data/data/Illumina_adapter_reads_PE_SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20 >> amtTrimmed.txt
done
```

Submitted batch job 1652405

After job is complete, make sure all data is there 

Check number of files 
```
ls -1 | wc -l 
```

64 files here

### 5) Rerun fastqc and multiqc

Count number of reads in trimmed files

```
nano read-counts_trim.sh
```

```
#!/bin/bash
#SBATCH -t 4:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@my.uri.edu
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Francois_data/scripts
#SBATCH --error="read-count_trim_error"
#SBATCH --output="read-count_trim"

for file in /data/putnamlab/jillashey/Francois_data/data/trimmed/*.trim.fq.gz
do
echo $file >> read-counts.txt
zcat $file | echo $((`wc -l`/4)) >> read-counts_trim.txt
done
```
Submitted batch job 1652530


Running fastqc on trimmed data

```
#!/bin/bash
#SBATCH -t 48:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=150GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@my.uri.edu
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Francois_data/scripts
#SBATCH --error="fastqc_out_trim_error"
#SBATCH --output="fastqc_out_trim"

module load FastQC/0.11.8-Java-1.8 #another version in Bluewaves if this one doesn't work

for file in /data/putnamlab/jillashey/Francois_data/data/trimmed/*.trim.fq.gz
do
fastqc $file --outdir /data/putnamlab/jillashey/Francois_data/fastqc_results/trimmed/
done

```

Submitted batch job 1652531

2x number of files so 128 fastqc files 

Run multiqc on trimmed fastqc files 

```
module load MultiQC/1.7-foss-2018b-Python-2.7.15
multiqc /data/putnamlab/jillashey/Francois_data/fastqc_results/trimmed/*trim_fastqc.zip
```

Secure copy multiqc output to local computer 

```
scp jillashey@bluewaves.uri.edu:/data/putnamlab/jillashey/Francois_data/multiqc/trimmed/multiqc_report.html /Users/jillashey/Desktop/
scp -r jillashey@bluewaves.uri.edu:/data/putnamlab/jillashey/Francois_data/multiqc/trimmed/multiqc_data /Users/jillashey/Desktop/MultiQC 
```

![fastqc_sequence_counts_plot.png](https://raw.githubusercontent.com/hputnam/Tufts_URI_CSM_RNASeq/master/_images/SenecaQC/trimmed/fastqc_sequence_counts_plot.png?token=APHKO3Y2WNWMQCEC2L5R3WK7EIHUS)

![fastqc_per_base_sequence_quality_plot.png](https://raw.githubusercontent.com/hputnam/Tufts_URI_CSM_RNASeq/master/_images/SenecaQC/trimmed/fastqc_per_base_sequence_quality_plot.png?token=APHKO35W2SK5X22OJCIWLDK7EIHYK)

![fastqc_per_sequence_quality_scores_plot.png](https://raw.githubusercontent.com/hputnam/Tufts_URI_CSM_RNASeq/master/_images/SenecaQC/trimmed/fastqc_per_sequence_quality_scores_plot.png?token=APHKO34U5XSYTXU6CVMNCRS7EIH3I)

![fastqc_per_base_sequence_content_plot.png](https://raw.githubusercontent.com/hputnam/Tufts_URI_CSM_RNASeq/master/_images/SenecaQC/trimmed/fastqc_per_base_sequence_content_plot.png?token=APHKO3ZIQ6GIESOAWT7ZQNK7EIIKE)

![fastqc_per_sequence_gc_content_plot.png](https://raw.githubusercontent.com/hputnam/Tufts_URI_CSM_RNASeq/master/_images/SenecaQC/trimmed/fastqc_per_sequence_gc_content_plot.png?token=APHKO34N6G7WSULPGCXJXBC7EIIPW)

![fastqc_sequence_duplication_levels_plot.png](https://raw.githubusercontent.com/hputnam/Tufts_URI_CSM_RNASeq/master/_images/SenecaQC/trimmed/fastqc_sequence_duplication_levels_plot.png?token=APHKO35WWX6VPHVXK6FYQL27EIISS)

![fastqc_overrepresented_sequencesi_plot.png](https://raw.githubusercontent.com/hputnam/Tufts_URI_CSM_RNASeq/master/_images/SenecaQC/trimmed/fastqc_overrepresented_sequencesi_plot.png?token=APHKO3522PRTQJG6GCP3L627EIIWA)

![fastqc_adapter_content_plot.png](https://raw.githubusercontent.com/hputnam/Tufts_URI_CSM_RNASeq/master/_images/SenecaQC/trimmed/fastqc_adapter_content_plot.png?token=APHKO326DYTYY4LOSNO2G3S7EIIYY)

### 6) Run STAR
```
#!/bin/bash
#SBATCH -t 12:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@my.uri.edu
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Francois_data/scripts
#SBATCH --error="gunzip_trim_out_error"
#SBATCH --output="gunzip_trim_out"

for file in /data/putnamlab/jillashey/Francois_data/data/trimmed/*.fq.gz
do
gunzip $file 
done
```

Submitted batch job 1652887

```
nano staralign.sh
```

```
#!/bin/bash
#SBATCH --mem=64G
module load STAR/2.5.3a-foss-2016b
STAR --runMode alignReads --quantMode TranscriptomeSAM --outTmpDir /data/putnamlab/jillashey/Francois_data/output/STAR/${SLURM_ARRAY_TASK_ID}_TMP --readFilesIn /data/putnamlab/jillashey/Francois_data/data/trimmed/${SLURM_ARRAY_TASK_ID}.fastq.trim.fq --genomeDir /data/putnamlab/jillashey/Connelly_data/output/STAR/GenomeIndex/ --sjdbGTFfeatureExon exon --sjdbGTFtagExonParentTranscript Parent --sjdbGTFfile /data/putnamlab/jillashey/Connelly_data/genome/GCF_003704095.1_ASM370409v1_genomic.gff --twopassMode Basic --twopass1readsN -1 --outStd Log BAM_Unsorted BAM_Quant --outSAMtype BAM Unsorted --outReadsUnmapped Fastx --outFileNamePrefix /data/putnamlab/jillashey/Francois_data/output/STAR/staralign_pdam/${SLURM_ARRAY_TASK_ID}.

# Submit as job 
sbatch --array staralign.sh
```