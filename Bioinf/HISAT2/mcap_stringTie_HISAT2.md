StringTie - hisat2, mcap genome

First, I was using Mcap.GFFannotation.fixed_all_exons_no_spaces.gff as the GFF file, but I may need to use a different one that has the correct ids (ie not all exon). Decided to use Mcap.GFFannotation.fixed_no_spaces.gff,
which is one that I fixed myself. I formatted column 9 (transcript_id=, gene_id=, etc), but kept names of ids

 as it is the original gff file and erin has used it succeessfully in the past 

a) Assemble reads with genome annotation

```
nano stringTie_mcap_assemble.sh

#!/bin/bash
#SBATCH --mem=64GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

# StringTie to assemble transcripts for each sample with the annotation file

F=/data/putnamlab/jillashey/Francois_data/stringTie_hisat2/mcap/BAMfiles/

array1=($(ls $F/*bam))
for i in ${array1[@]}; do
	stringtie -G /data/putnamlab/jillashey/genome/Mcap/Mcap.GFFannotation.fixed_no_spaces.gff -e -o ${i}.gtf ${i}
	echo "${i}"
done

sbatch stringTie_mcap_assemble.sh
```
Submitted batch job 1671941


b) Merge stringTie gtf results 

```
mv *gtf ../GTFfiles/

ls *gtf > mcap_mergelist.txt
cat mcap_mergelist.txt 

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

stringtie --merge -p 8 -G /data/putnamlab/jillashey/genome/Mcap/Mcap.GFFannotation.fixed_no_spaces.gff -o stringtie_mcap_merged.gtf mcap_mergelist.txt

```

c) Assess assembly quality

```
module load gffcompare/0.11.5-foss-2018b

gffcompare -r /data/putnamlab/jillashey/genome/Mcap/Mcap.GFFannotation.fixed_no_spaces.gff -o ../merged stringtie_mcap_merged.gtf

9352 reference transcripts loaded.
  9352 query transfrags loaded.
```

This number seems very low and its weird that they match. when checking merged.stats, sensitivity and precision are all at 100%, which isn't possible. Not sure if its the gff file or my code?

I'm going to continue, but unsure about results

d) Re-estimate assembly

```
nano stringTie_mcap_re-assemble.sh

#!/bin/bash
#SBATCH --mem=64GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/stringTie_hisat2/mcap/BAMfiles/

array1=($(ls $F/*bam))
for i in ${array1[@]}
do
stringtie -e -G /data/putnamlab/jillashey/genome/Mcap/Mcap.GFFannotation.fixed_no_spaces.gff -o ${i}.merge.gtf ${i}
echo "${i}"
done

sbatch stringTie_mcap_re-assemble.sh
```

e) Create gene matrix

```
module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/stringTie_hisat2/mcap/GTFfiles_merged/

array2=($(ls *merge.gtf))

for i in ${array2[@]}
do
echo "${i} $F${i}" >> sample_list_mcap.txt
done

python prepDE.py -g gene_count_mcap_matrix.csv -i sample_list_mcap.txt

```
