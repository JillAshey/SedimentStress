### StringTie 

STAR, pdam genome - Using annotation file (GFFS) and fasta file (scaffold) from http://pdam.reefgenomics.org/download/ 


a) Assemble reads with genome annotation

```
nano stringTie_pdam_rgGFF_assemble.sh

#!/bin/bash
#SBATCH --mem=64GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

# StringTie to assemble transcripts for each sample with the annotation file

F=/data/putnamlab/jillashey/Francois_data/stringTie_star/pdam_rgGFF

array1=($(ls $F/*bam))
for i in ${array1[@]}; do
	stringtie -G /data/putnamlab/jillashey/genome/Pdam/pdam_annotation_AddTranscript_id_fixed.gtf -e -o ${i}.gtf ${i}
	echo "${i}"
done

sbatch stringTie_pdam_rgGFF_assemble.sh
```
Submitted batch job 1672895

b) Merge stringTie gtf results 

```
# Move only pdam samples
mv *gtf ../GTFfiles

ls *gtf > pdam_rgGFF_mergelist.txt
cat pdam_rgGFF_mergelist.txt 

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

stringtie --merge -p 8 -G /data/putnamlab/jillashey/genome/Pdam/pdam_annotation_AddTranscript_id_fixed.gtf -o stringtie_pdam_rgGFF_merged.gtf pdam_rgGFF_mergelist.txt

```

Gives me this warning: not sure why, but this was also given to me when 
```
Warning: invalid start coordinate at line:
###			NA	NA				
Warning: invalid start coordinate at line:
###			NA	NA
```

c) Assess assembly quality

```
module load gffcompare/0.11.5-foss-2018b

gffcompare -r /data/putnamlab/jillashey/genome/Pdam/pdam_annotation_AddTranscript_id_fixed.gtf -o merged stringtie_pdam_rgGFF_merged.gtf

1675 reference transcripts loaded.
  1675 query transfrags loaded.

```

merge stats say everything is 100% ? I know its supposed to be high, but is it supposed to be perfect like that?

d) Re-estimate assembly

```
nano stringTie_pdam_rgGFF_re-assemble.sh

#!/bin/bash
#SBATCH --mem=64GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/stringTie_hisat2/pdam/BAMfiles

array1=($(ls $F/*bam))
for i in ${array1[@]}
do
stringtie -e -G /data/putnamlab/jillashey/genome/Pdam/pdam_annotation_AddTranscript_id_fixed.gtf -o ${i}.merge.gtf ${i}
echo "${i}"
done

sbatch stringTie_pdam_rgGFF_re-assemble.sh
```

Submitted batch job 1672896

e) Create gene matrix

```
module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/stringTie_star/pdam_rgGFF/GTFfiles_merged/

array2=($(ls *merge.gtf))

for i in ${array2[@]}
do
echo "${i} $F${i}" >> sample_list_pdam_rgGFF.txt
done

python prepDE.py -g gene_count_pdam_rgGFF_matrix.csv -i sample_list_pdam_rgGFF.txt

```

f) Secure-copy gene counts onto local computer 

```
scp jillashey@bluewaves.uri.edu:/data/putnamlab/jillashey/Francois_data/stringTie_star/pdam_rgGFF/GTFfiles_merged/gene_count_pdam_rgGFF_matrix.csv /Users/jillashey/Desktop/
```
