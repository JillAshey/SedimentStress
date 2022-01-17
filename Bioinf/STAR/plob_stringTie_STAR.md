StringTie - STAR, plobata samples with plutea genome 

STAR was run without GFF file included 

a) Move reads (aligned by Coordinate) from output to stringTie folder 

For some reason, 38_2 not running. That's okay for now, as its a Pdam sample.

```
mv *Aligned.sortedByCoord.out.bam ../../../stringTie_star/plob/
```

b) Assemble reads with genome annotation

```
nano stringTie_plob_assemble.sh

#!/bin/bash
#SBATCH --mem=64GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --error="assemble_plob_out_error"
#SBATCH --output="assemble_plob_out"

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

# StringTie to assemble transcripts for each sample with the annotation file

F=/data/putnamlab/jillashey/Francois_data/Hawaii/stringTie_star/plob/BAM

array1=($(ls $F/*bam))
for i in ${array1[@]}; do
	stringtie -G /data/putnamlab/jillashey/genome/Plutea/Plut.GFFannotation.fixed_transcript.gff -e -o ${i}.gtf ${i}
	echo "${i}"
done

sbatch stringTie_plob_assemble.sh
```

Submitted batch job 1768981

c) Merge stringTie gtf results 

At this point, I'm only going to move forward with Plob samples

```
# move only plob samples
mv 6_1.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.gtf 6_2.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.gtf 7_1.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.gtf 8_1.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.gtf 9_1.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.gtf 9_2.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.gtf 21_1.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.gtf 22_1.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.gtf 23_1.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.gtf 25_1.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.gtf 25_2.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.gtf 26_1.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.gtf 26_2.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.gtf 27_1.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.gtf 29_1.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.gtf 34_1.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.gtf ../GTF 
```

```

ls *gtf > plob_mergelist.txt
cat plob_mergelist.txt

module load StringTie/2.1.1-GCCcore-7.3.0

stringtie --merge -p 8 -G /data/putnamlab/jillashey/genome/Plutea/Plut.GFFannotation.fixed_transcript.gff -o stringtie_plob_merged.gtf plob_mergelist.txt
```

d) Assess assembly quality

```
module load gffcompare/0.11.5-foss-2018b

gffcompare -r /data/putnamlab/jillashey/genome/Plutea/Plut.GFFannotation.fixed_transcript.gff -o Plob.merged stringtie_plob_merged.gtf

  31126 reference transcripts loaded.
  31126 query transfrags loaded.
```

e) Re-estimate assembly 

```
nano stringTie_plob_re-assemble.sh

#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="re-assemble_plob_out_error"
#SBATCH --output="re-assemble_plob_out"

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/Hawaii/stringTie_star/plob/BAM # do i need a / after BAM?

array1=($(ls $F/*bam))
for i in ${array1[@]}
do
stringtie -e -G /data/putnamlab/jillashey/genome/Plutea/Plut.GFFannotation.fixed_transcript.gff -o ${i}.merge.gtf ${i}
echo "${i}"
done

sbatch stringTie_plob_re-assemble.sh
```

Submitted batch job 1768989

```
# move only plob samples
mv 6_1.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.merge.gtf 6_2.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.merge.gtf 7_1.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.merge.gtf 8_1.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.merge.gtf 9_1.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.merge.gtf 9_2.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.merge.gtf 21_1.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.merge.gtf 22_1.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.merge.gtf 23_1.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.merge.gtf 25_1.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.merge.gtf 25_2.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.merge.gtf 26_1.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.merge.gtf 26_2.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.merge.gtf 27_1.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.merge.gtf 29_1.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.merge.gtf 34_1.fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.merge.gtf ../GTF_merge 
```

f) Create gene matrix

```
module load StringTie/2.1.1-GCCcore-7.3.0
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/Hawaii/stringTie_star/plob/GTF_merge/

array2=($(ls *merge.gtf))

for i in ${array2[@]}
do
echo "${i} $F${i}" >> sample_list_plob.txt
done

python prepDE.py -g gene_count_plob_only_matrix.csv -i sample_list_plob.txt
```

g) Secure-copy gene counts onto local computer

```
scp jillashey@bluewaves.uri.edu:/data/putnamlab/jillashey/Francois_data/Hawaii/stringTie_star/plob/GTF_merge/gene_count_plob_only_matrix.csv /Users/jillashey/Desktop/Putnamlab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/
```

20220116

Going to try to run re-assemble step with the stringtie_plob_merged.gtf file instead of the Plut.GFFannotation.fixed_transcript.gff file. 

```
cd /data/putnamlab/jillashey/Francois_data/Hawaii/stringTie_star/plob/BAM/plob_only

nano stringTie_plob_re-assemble_merge.sh

#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="re-assemble_merge_plob_out_error"
#SBATCH --output="re-assemble_merge_plob_out"

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/Hawaii/stringTie_star/plob/BAM/plob_only # do i need a / after BAM?

array1=($(ls $F/*bam))
for i in ${array1[@]}
do
stringtie -e -G /data/putnamlab/jillashey/Francois_data/Hawaii/stringTie_star/plob/GTF/stringtie_plob_merged.gtf -o ${i}.merge.gtf ${i}
echo "${i}"
done

sbatch stringTie_plob_re-assemble_merge.sh
```

Submitted batch job 1960767

Move samples to new GTF folder 

f) Create gene matrix

```
module load StringTie/2.1.1-GCCcore-7.3.0
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/Hawaii/stringTie_star/plob/GTF_merge_20220116/

array2=($(ls *merge.gtf))

for i in ${array2[@]}
do
echo "${i} $F${i}" >> sample_list_plob.txt
done

python prepDE.py -g gene_count_plob_matrix_merge.csv -i sample_list_plob.txt

```

I'm sticking w/ the original gene count matrix that I made 