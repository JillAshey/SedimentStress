StringTie - hisat2, pcomp genome

Not sure which gff annotation file is best to use

a) Assemble reads with genome annotation

```
nano stringTie_pcomp_assemble.sh

#!/bin/bash
#SBATCH --mem=64GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

# StringTie to assemble transcripts for each sample with the annotation file

F=/data/putnamlab/jillashey/Francois_data/stringTie_hisat2/pcomp/BAMfiles/

array1=($(ls $F/*bam))
for i in ${array1[@]}; do
	stringtie -G /data/putnamlab/jillashey/genome/Pcomp/Pcomp.GFFannotation.fixed.gff -e -o ${i}.gtf ${i}
	echo "${i}"
done

sbatch stringTie_pcomp_assemble.sh
```

Submitted batch job 1671954

b) Merge stringTie gtf results 

```
mv *gtf ../GTFfiles/

ls *gtf > pcomp_mergelist.txt
cat pcomp_mergelist.txt 

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

stringtie --merge -p 8 -G /data/putnamlab/jillashey/genome/Pcomp/Pcomp.GFFannotation.fixed.gff -o stringtie_pcomp_merged.gtf pcomp_mergelist.txt

```

c) Assess assembly quality

```
module load gffcompare/0.11.5-foss-2018b

gffcompare -r /data/putnamlab/jillashey/genome/Pcomp/Pcomp.GFFannotation.fixed.gff -o ../merged stringtie_pcomp_merged.gtf

74728 reference transcripts loaded.
  11 duplicate reference transcripts discarded.
  74727 query transfrags loaded.
```
Merged stats output has sensitivity and precision close to 100% for all levels, which seems impossible. may be the annotation file that is creating this weird stuff 

d) Re-estimate assembly 

```
nano stringTie_pcomp_re-assemble.sh

#!/bin/bash
#SBATCH --mem=64GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/stringTie_hisat2/pcomp/BAMfiles

array1=($(ls $F/*bam))
for i in ${array1[@]}
do
stringtie -e -G /data/putnamlab/jillashey/genome/Pcomp/Pcomp.GFFannotation.fixed.gff -o ${i}.merge.gtf ${i}
echo "${i}"
done

sbatch stringTie_pcomp_re-assemble.sh
```

Submitted batch job 1671960

e) Create gene matrix

```
module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/stringTie_hisat2/pcomp/GTFfiles_merge/

array2=($(ls *merge.gtf))

for i in ${array2[@]}
do
echo "${i} $F${i}" >> sample_list_pcomp.txt
done

python prepDE.py -g gene_count_pcomp_matrix.csv -i sample_list_pcomp.txt

```
