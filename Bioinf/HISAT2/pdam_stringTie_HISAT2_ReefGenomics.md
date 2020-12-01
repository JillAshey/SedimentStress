### StringTie 

HISAT2, pdam genome - Using annotation file (GFFS) and fasta file (scaffold) from http://pdam.reefgenomics.org/download/ 

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

F=/data/putnamlab/jillashey/Francois_data/stringTie_hisat2/pdam_rgGFF

array1=($(ls $F/*bam))
for i in ${array1[@]}; do
	stringtie -G /data/putnamlab/jillashey/genome/Pdam/pdam_annotation_AddTranscript_id_fixed.gtf -e -o ${i}.gtf ${i}
	echo "${i}"
done

sbatch stringTie_pdam_rgGFF_assemble.sh
```
Submitted batch job 1672986

b) Merge stringTie gtf results 

```
# Move samples
mv *gtf ../GTFfiles

ls *gtf > pdam_rgGFF_hisat2_mergelist.txt
cat pdam_rgGFF_hisat2_mergelist.txt 

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

stringtie --merge -p 8 -G /data/putnamlab/jillashey/genome/Pdam/pdam_annotation_AddTranscript_id_fixed.gtf -o stringtie_pdam_rgGFF_hisat2_merged.gtf pdam_rgGFF_hisat2_mergelist.txt

```

Again got this error: 

```
Warning: invalid start coordinate at line:
###			NA	NA				
Warning: invalid start coordinate at line:
###			NA	NA
```
Not sure why???? I think it must be something with the format of the reef genomics annotation file. Doesn't seem to be affected anything in the analysis at first glance, but will keep it in mind as I move through downstream analysis 

c) Assess assembly quality

```
module load gffcompare/0.11.5-foss-2018b

gffcompare -r /data/putnamlab/jillashey/genome/Pdam/pdam_annotation_AddTranscript_id_fixed.gtf -o merged stringtie_pdam_rgGFF_hisat2_merged.gtf

1675 reference transcripts loaded.
  1675 query transfrags loaded.
```

Interesting, same results loaded as pdam with STAR with the reef genomics annotation file. Is this due to the annotation file? 

d) Re-estimate assembly

```
nano stringTie_pdam_rgGFF_hisat2_re-assemble.sh

#!/bin/bash
#SBATCH --mem=64GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/stringTie_hisat2/pdam_rgGFF/BAMfiles

array1=($(ls $F/*bam))
for i in ${array1[@]}
do
stringtie -e -G /data/putnamlab/jillashey/genome/Pdam/pdam_annotation_AddTranscript_id_fixed.gtf -o ${i}.merge.gtf ${i}
echo "${i}"
done

sbatch stringTie_pdam_rgGFF_hisat2_re-assemble.sh
```

Submitted batch job 1673013

e) Create gene matrix

```
module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/stringTie_hisat2/pdam_rgGFF/GTFfiles_merged/

array2=($(ls *merge.gtf))

for i in ${array2[@]}
do
echo "${i} $F${i}" >> sample_list_pdam_rgGFF_hisat2.txt
done

python prepDE.py -g gene_count_pdam_rgGFF_hisat2_matrix.csv -i sample_list_pdam_rgGFF_hisat2.txt
```

Xxxxxxxxxxxxx
### Using original gff3 file from RG
In the code above, I used the gff file that I modified when running star. I added transcript_id= to the last column. Having some issues comparing to Polina's output and she used the original gff file, so I am going to use that one below

a) Assemble reads with genome annotation

```
nano stringTie_pdam_rgGFF_original_assemble.sh

#!/bin/bash
#SBATCH --mem=64GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

# StringTie to assemble transcripts for each sample with the annotation file

F=/data/putnamlab/jillashey/Francois_data/stringTie_hisat2/pdam_rgGFF/BAMfiles

array1=($(ls $F/*bam))
for i in ${array1[@]}; do
	stringtie -G /data/putnamlab/jillashey/genome/Pdam/pdam_annotation.gff3 -e -o ${i}.OG.gtf ${i}
	echo "${i}"
done

sbatch stringTie_pdam_rgGFF_original_assemble.sh
```
Submitted batch job 1673884

b) Merge stringTie gtf results 

```
# Move samples
mv *gtf ../GTFfiles/originalGFF

ls *gtf > pdam_rgGFF_original_hisat2_mergelist.txt
cat pdam_rgGFF_original_hisat2_mergelist.txt 

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

stringtie --merge -p 8 -G /data/putnamlab/jillashey/genome/Pdam/pdam_annotation.gff3 -o stringtie_pdam_original_rgGFF_hisat2_merged.gtf pdam_rgGFF_original_hisat2_mergelist.txt

```

c) Assess assembly quality

```
module load gffcompare/0.11.5-foss-2018b

gffcompare -r /data/putnamlab/jillashey/genome/Pdam/pdam_annotation.gff3 -o original_merged stringtie_pdam_original_rgGFF_hisat2_merged.gtf

26077 reference transcripts loaded.
  26077 query transfrags loaded.
# looks good
```

Not getting the warning errors that I got above which is good

d) Re-estimate assembly

```
nano stringTie_pdam_rgGFF_original_hisat2_re-assemble.sh

#!/bin/bash
#SBATCH --mem=64GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/stringTie_hisat2/pdam_rgGFF/BAMfiles

array1=($(ls $F/*bam))
for i in ${array1[@]}
do
stringtie -e -G /data/putnamlab/jillashey/genome/Pdam/pdam_annotation.gff3 -o ${i}.original_merge.gtf ${i}
echo "${i}"
done

sbatch stringTie_pdam_rgGFF_original_hisat2_re-assemble.sh
```

Submitted batch job 1673920

e) Create gene matrix

```
module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/stringTie_hisat2/pdam_rgGFF/GTFfiles_merged/original_gff/

array2=($(ls *merge.gtf))

for i in ${array2[@]}
do
echo "${i} $F${i}" >> sample_list_pdam_rgGFF_original_hisat2.txt
done

python prepDE.py -g gene_count_pdam_rgGFF_original_hisat2_matrix.csv -i sample_list_pdam_rgGFF_original_hisat2.txt

# gives me this error
Error: could not locate transcript pdam_00004790-RA entry for sample 11_1.fastq.trim.fq.bam.original_merge.gtf
Traceback (most recent call last):
  File "prepDE.py", line 281, in <module>
    geneDict.setdefault(geneIDs[i],{}) #gene_id
KeyError: 'pdam_00004790-RA'
```