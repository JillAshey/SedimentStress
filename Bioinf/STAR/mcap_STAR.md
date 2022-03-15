## STAR M. capitata

Date: 20220217

For sediment stress paper, we initially did not have Mcap metadata and so did not analyze Mcap samples. Now metadata is found so I need to re-run STAR using the Mcap version2 genomic info 

Deleted old Mcap genome directory off of server bc those files were all of old Mcap version. Now using files from http://cyanophora.rutgers.edu/montipora/ (version 2). Download and unzip all files before moving forward.

File quality and trimming has already been done for these reads, so will go right into STAR

a) After downloading onto server, unzip files 

```
ls
Montipora_capitata_HIv2.assembly.fasta.gz  Montipora_capitata_HIv2.genes.cds.fna.gz  Montipora_capitata_HIv2.genes.gff3.gz  Montipora_capitata_HIv2.genes.pep.faa.gz

gunzip *
```

Count number of seqs per file 

```
# Protein
zgrep -c ">" Montipora_capitata_HIv2.genes.pep.faa 
55086

# Transcript
zgrep -c ">" Montipora_capitata_HIv2.genes.cds.fna 
55086
```

b) Generate genome index

STAR needs transcript_id identifier to run properly, so added that identifier for Pacuta in R. Renamed gff file as Mcap.gff.annotations.fixed_transcript.gff3 and uploaded to server.

Moved all old STAR output into new directory called 


```
cd /data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR
mkdir Mcap_V1
mv *mcap* Mcap_V1/
```

Run genome index step of STAR

Write as a script and submit as job 

```
mkdir GenomeIndex_mcapV2

nano star_index_mcap.sh

#!/bin/bash
#SBATCH -t 48:00:00
#SBATCH --nodes=1 --ntasks-per-node=2
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Francois_data/Hawaii/scripts
#SBATCH --error="star_index_mcap_out_error"
#SBATCH --output="star_index_mcap_out"

module load STAR/2.5.3a-foss-2016b # on bluewaves

STAR --runThreadN 10 --runMode genomeGenerate --genomeDir /data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/GenomeIndex_mcapV2 --genomeFastaFiles /data/putnamlab/jillashey/genome/Mcap/Montipora_capitata_HIv2.assembly.fasta --sjdbGTFfile /data/putnamlab/jillashey/genome/Mcap/Mcap.gff.annotations.fixed_transcript.gff3

sbatch star_index_mcap.sh
```

Submitted batch job 1978534

20220219 - okay so it ran, but something got messed up. there is no info in the genome index for sjdbList files. It might be because of how I separated out the transcript_id in R? Going to go back now and check code

20220221 - I separated out the transcript_id in R using a different method, copied the new gff fixed file to server, and now going to re-run ```star_index_mcap.sh```

Submitted batch job 1979176

It seems to have worked! Good amount of info in the sjdbList files. Now to align reads against mcap

c) Align reads 

The mcap samples are: `5_2.fastq.gz, 10_2.fastq.gz, 13_2.fastq.gz, 14_2.fastq.gz, 15_2.fastq.gz, 16_2.fastq.gz, 19_2.fastq.gz, 20_2.fastq.gz, 40_2.fastq.gz, 45_2.fastq.gz, 48_2.fastq.gz`. Going to separate them out to run on STAR first and then will run all other samples


Link to trimmed files 

```
mkdir AlignReads_mcapV2
cd AlignReads_mcapV2
ln -s /data/putnamlab/jillashey/Francois_data/Hawaii/data/trimmed/*trim.fq .
```

Move mcap samples to their own directory 

```
mkdir mcap_only
mv 5_2.fastq.trim.fq 10_2.fastq.trim.fq 13_2.fastq.trim.fq 14_2.fastq.trim.fq 15_2.fastq.trim.fq 16_2.fastq.trim.fq 19_2.fastq.trim.fq 20_2.fastq.trim.fq 40_2.fastq.trim.fq 45_2.fastq.trim.fq 48_2.fastq.trim.fq mcap_only/
cd mcap_only
```

Let's do a test sample first 

```
nano test_5_2_align.sh

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

STAR --runMode alignReads --quantMode TranscriptomeSAM --outTmpDir /data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/test_TMP --readFilesIn 5_2.fastq.trim.fq --genomeDir /data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/GenomeIndex_mcapV2/ --twopassMode Basic --twopass1readsN -1 --outStd Log BAM_Unsorted BAM_Quant --outSAMtype BAM Unsorted SortedByCoordinate --outReadsUnmapped Fastx --outFileNamePrefix 5_2_test.

sbatch test_5_2_align.sh
```

Submitted batch job 1979184

Keep getting this error: EXITING because of fatal ERROR: could not make temporary directory: /data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/
SOLUTION: (i) please check the path and writing permissions 
 (ii) if you specified --outTmpDir, and this directory exists - please remove it before running STAR

Feb 21 23:43:57 ...... FATAL ERROR, exiting

Not sure why... going to try later 

Okay tried again. I had the --outTmpDir going to `/data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/AlignReads_mcapV2` but changed it so it would go to `/data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/test_TMP `. Submitted batch job 1979524. Success! Took ~75 mins 

Now will run w/ all mcap samples 

```
nano AlignReads_mcap_only.sh

#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="Align_mcap_only_out_error"
#SBATCH --output="Align_mcap_only_out"

module load STAR/2.5.3a-foss-2016b

F=/data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/AlignReads_mcapV2/mcap_only

array1=($(ls $F/*trim.fq))
for i in ${array1[@]}
do
STAR --runMode alignReads --quantMode TranscriptomeSAM --outTmpDir ${i}_TMP --readFilesIn ${i} --genomeDir /data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/GenomeIndex_mcapV2/ --twopassMode Basic --twopass1readsN -1 --outStd Log BAM_Unsorted BAM_Quant --outSAMtype BAM Unsorted SortedByCoordinate --outReadsUnmapped Fastx --outFileNamePrefix ${i}.
done 

sbatch --array=1-11 AlignReads_mcap_only.sh # submit as an array job 
```

Submitted batch job 1979526

Cool it ran in a few hours. 3 samples still running. Now will align the rest of samples against mcap genome 

```
nano AlignReads_mcap.sh

#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="Align_mcap_out_error"
#SBATCH --output="Align_only_out"

module load STAR/2.5.3a-foss-2016b

F=/data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/AlignReads_mcapV2

array1=($(ls $F/*trim.fq))
for i in ${array1[@]}
do
STAR --runMode alignReads --quantMode TranscriptomeSAM --outTmpDir ${i}_TMP --readFilesIn ${i} --genomeDir /data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/GenomeIndex_mcapV2/ --twopassMode Basic --twopass1readsN -1 --outStd Log BAM_Unsorted BAM_Quant --outSAMtype BAM Unsorted SortedByCoordinate --outReadsUnmapped Fastx --outFileNamePrefix ${i}.
done 

sbatch --array=1-53%15 AlignReads_mcap.sh # submit as an array job 
```

Submitted batch job 1979560

d) Perform gene counts with stringTie

Make folders for stringtie results 

```
cd /data/putnamlab/jillashey/Francois_data/Hawaii/stringTie_star
mkdir mcap_V2
cd mcap_V2
mkdir BAM GTF GTF_merge
```

Move BAM files to stringTie folder 

```
cd /data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/AlignReads_mcapV2/mcap_only
mv *Aligned.sortedByCoord.out.bam /data/putnamlab/jillashey/Francois_data/Hawaii/stringTie_star/mcap_V2/BAM
```

For some reason the 14_2, 19_2, and 45_2 BAM files is empty? so something must have gone wrong in STAR, but it still gave me % mapping and all that. going back to STAR to rerun. 

When checking file sizes with `ls -l`, it looks like all the files in the STAR folder are okay except for 45_2. So why is 14_2 and 19_2 empty? Weird. it looks like some samples just take longer to align in STAR? Maybe because that particular sample has a lot of reads, but looking at raw read counts, thats not necessarily true. Maybe STAR just quit early? Maybe its an issue with running jobs as an array. idk. this is not ideal, but im just going to align them with their own scripts 

14_2 alignment

```
nano 14_2_align.sh

#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="14_2_align_mcap_out_error"
#SBATCH --output="14_2_align_mcap_out_error"

module load STAR/2.5.3a-foss-2016b

STAR --runMode alignReads --quantMode TranscriptomeSAM --outTmpDir /data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/test_TMP --readFilesIn 14_2.fastq.trim.fq --genomeDir /data/putnamlab/jillashey/Francois_data/Hawaii/output/STAR/GenomeIndex_mcapV2/ --twopassMode Basic --twopass1readsN -1 --outStd Log BAM_Unsorted BAM_Quant --outSAMtype BAM Unsorted SortedByCoordinate --outReadsUnmapped Fastx --outFileNamePrefix 14_2.

sbatch 14_2_align.sh
```

Submitted batch job 1979568 - canceled, was still pending. waiting for other STAR alignment to finish running 

also may need to do 48_2? Very confusing...when I look at the Log.final.out file, it gives me % mapped and all the other stats. But when I look at the Log.progress.out file, it doesn't finish mapping...But when I look at the Log.out file, it says all done at the bottom.

Samples that did not finish mapping (according to Log.progress.out file): 15_2, 20_2,

40_2 doesnt even have any data in Log.progress.out file. Ugh in the Log.out file, it only got to loading the genome. How does it have a Log.final.out file then?

okay just let everyhing run. looks like everything was just still running and i was jumping the gun moving forward w/ the analysis. 

NOW assemble and estimate reads 

Assemble reads w/ genome annotation

```
cd /data/putnamlab/jillashey/Francois_data/Hawaii/stringTie_star/mcap_V2/BAM

nano stringTie_mcap_assemble.sh

#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=2
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="assemble_mcap_out_error"
#SBATCH --output="assemble_mcap_out"

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/Hawaii/stringTie_star/mcap_V2/BAM

array1=($(ls $F/*bam))
for i in ${array1[@]}; do
	stringtie -G /data/putnamlab/jillashey/genome/Mcap/Mcap.gff.annotations.fixed_transcript.gff3 -e -o ${i}.gtf ${i}
	echo "${i}"
done

sbatch stringTie_mcap_assemble.sh
```

Submitted batch job 1981735. Finished in ~1hr

Merge stringTie gtf results 

```
# move gtf files to their own folder 
cd /data/putnamlab/jillashey/Francois_data/Hawaii/stringTie_star/mcap_V2/BAM
mv *gtf ../GTF
cd ../GTF

ls *gtf > mcap_mergelist.txt
cat mcap_mergelist.txt
module load StringTie/2.1.1-GCCcore-7.3.0
stringtie --merge -p 8 -G /data/putnamlab/jillashey/genome/Mcap/Mcap.gff.annotations.fixed_transcript.gff3 -o stringtie_mcap_merged.gtf mcap_mergelist.txt
```

Assess assembly quality 

```
module load gffcompare/0.11.5-foss-2018b

gffcompare -r /data/putnamlab/jillashey/genome/Mcap/Mcap.gff.annotations.fixed_transcript.gff3 -o Mcap.merged stringtie_mcap_merged.gtf
	55086 reference transcripts loaded.
	55086 query transfrags loaded.
# good, same # of protein seqs
```

Re-estimate assembly 

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

F=/data/putnamlab/jillashey/Francois_data/Hawaii/stringTie_star/mcap_V2/BAM

array1=($(ls $F/*bam))
for i in ${array1[@]}
do
stringtie -e -G /data/putnamlab/jillashey/genome/Mcap/Mcap.gff.annotations.fixed_transcript.gff3 -o ${i}.merge.gtf ${i}
echo "${i}"
done

sbatch stringTie_mcap_re-assemble.sh
```

Submitted batch job 1981736

Create gene matrix 

```
# copy python script to mcap folder 
cd /data/putnamlab/jillashey/Francois_data/Hawaii/stringTie_star/pacuta/GTF_merge
cp prepDE.py ../../mcap_V2/GTF_merge/

# move merge.gtf files to their own folder 
cd /data/putnamlab/jillashey/Francois_data/Hawaii/stringTie_star/mcap_V2/BAM
mv * merge.gtf ../GTF_merge
cd ../GTF_merge

# make gene matrix 
module load StringTie/2.1.1-GCCcore-7.3.0
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/Hawaii/stringTie_star/mcap_V2/GTF_merge/

array2=($(ls *merge.gtf))

for i in ${array2[@]}
do
echo "${i} $F${i}" >> sample_list_mcap.txt
done

python prepDE.py -g gene_count_mcap_matrix.csv -i sample_list_mcap.txt
```

HOORAY!!!!!!! FINALLY a gene count matrix for mcap!!!

Secure-copy gene counts onto local computer

```
scp jillashey@bluewaves.uri.edu:/data/putnamlab/jillashey/Francois_data/Hawaii/stringTie_star/mcap_V2/GTF_merge/gene_count_mcap_matrix.csv /Users/jillashey/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/mcap
```