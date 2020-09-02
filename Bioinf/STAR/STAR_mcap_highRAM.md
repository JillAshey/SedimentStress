### STAR - high RAM time

Talked with Kate today, super helpful. she recommended that I do a bunch of different things: 

1) Use newest version of STAR - need to email Kevin Bryan about getting newest version of STAR (2.7.5c)

2) give STAR all the memory I can (```--mem=510G```) - more RAM

3) Check if genome is unusually small or fragmented. If it is, there are arguments in STAR to deal with that

4) Find new Pcomp genome & info. Mine looks corrupt and the new one from the website doesn't have all the info I need



#### Increased memory 

STAR is RAM-intensive, which might explain why my samples are not being processed all the way. I'm going to run a sample with a lot of memory and see what happens. 

a) Generate genome index

```
mkdir GenomeIndex_mcap_highRAM

module load STAR/2.5.3a-foss-2016b

STAR --runThreadN 10 --runMode genomeGenerate --genomeDir /data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/GenomeIndex_mcap_highRAM --genomeFastaFiles /data/putnamlab/jillashey/genome/Mcap/Mcap.genome_assembly.fa --sjdbGTFfile /data/putnamlab/jillashey/genome/Mcap/Mcap.GFFannotation.fixed_all_exons_no_spaces.gff

```

b) Align reads 

I'm going to start aligning just one sample with a lot of RAM

40_2

```
mkdir AlignReads_mcap_highRAM
cd AlignReads_mcap_highRAM

nano AlignReads_mcap_40_2_highRAM.sh

#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=250G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="Align_mcap_40_2_out_error"
#SBATCH --output="Align_mcap_40_2_out"

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads \
--quantMode TranscriptomeSAM \
--outTmpDir /data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/AlignReads_mcap_highRAM/40_2.fastq.trim.fq_TMP \
--readFilesIn /data/putnamlab/jillashey/Francois_data/Hawaii/data/trimmed/40_2.fastq.trim.fq \
--genomeDir /data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/GenomeIndex_mcap_highRAM/ \
--twopassMode Basic \
--twopass1readsN -1 \
--outStd Log BAM_Unsorted BAM_Quant \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFileNamePrefix /data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/AlignReads_mcap_highRAM/40_2.fastq.trim.fq.

sbatch AlignReads_mcap_40_2_highRAM.sh

```

Submitted batch job 1692785 -- may take a while because needs lots of resources to run
