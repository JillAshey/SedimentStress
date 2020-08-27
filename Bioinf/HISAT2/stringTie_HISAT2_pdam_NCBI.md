### StringTie - hisat2, pdam genome from NCBI


a) Assemble reads with genome annotation

```
nano stringTie_pdam_assemble.sh

#!/bin/bash
#SBATCH --mem=64GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

# StringTie to assemble transcripts for each sample with the annotation file

F=/data/putnamlab/jillashey/Francois_data/stringTie_hisat2/pdam/BAMfiles/

array1=($(ls $F/*bam))
for i in ${array1[@]}; do
	stringtie -G /data/putnamlab/jillashey/genome/Pdam/GCF_003704095.1_ASM370409v1_genomic.gff -e -o ${i}.gtf ${i}
	echo "${i}"
done

sbatch stringTie_pdam_assemble.sh
```

Submitted batch job 1673031

b) Merge stringTie gtf results 

```
# Move only pdam samples
mv *gtf ../GTFfiles/

ls *gtf > pdam_mergelist.txt
cat pdam_mergelist.txt 

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

stringtie --merge -p 8 -G /data/putnamlab/jillashey/genome/Pdam/GCF_003704095.1_ASM370409v1_genomic.gff -o stringtie_pdam_merged.gtf pdam_mergelist.txt

```

c) Assess assembly quality

```
module load gffcompare/0.11.5-foss-2018b

gffcompare -r /data/putnamlab/jillashey/genome/Pdam/GCF_003704095.1_ASM370409v1_genomic.gff -o ../merged stringtie_pdam_merged.gtf

28862 reference transcripts loaded.
  28862 query transfrags loaded.

```

d) Re-estimate assembly

```
nano stringTie_pdam_re-assemble.sh

#!/bin/bash
#SBATCH --mem=64GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/stringTie_hisat2/pdam/BAMfiles/

array1=($(ls $F/*bam))
for i in ${array1[@]}
do
stringtie -e -G /data/putnamlab/jillashey/genome/Pdam/GCF_003704095.1_ASM370409v1_genomic.gff -o ${i}.merge.gtf ${i}
echo "${i}"
done

sbatch stringTie_pdam_re-assemble.sh
```

Submitted batch job 1673045

e) Create gene matrix

```
module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/stringTie_hisat2/pdam/GTFfiles_merged

array2=($(ls *merge.gtf))

for i in ${array2[@]}
do
echo "${i} $F/${i}" >> sample_list_pdam.txt
done

python prepDE.py -g gene_count_pdam_NCBI_hisat2_matrix.csv -i sample_list_pdam.txt

```

f) Secure-copy gene counts onto local computer 

```
scp jillashey@bluewaves.uri.edu:/data/putnamlab/jillashey/Francois_data/stringTie_hisat2/pdam/GTFfiles_merged/gene_count_pdam_NCBI_hisat2_matrix.csv /Users/jillashey/Desktop/
```









### Another gff annotation file used below 

I found that the Tufts server has a slightly different GFF file than I have (and that is on the NCBI server). The one of the Tufts server is better annotated and was done with maker. Not sure where it came from.


a) Assemble reads with genome annotation

```
nano stringTie_pdamSamplesOnly_assemble_TuftsGFF.sh

#!/bin/bash
#SBATCH --mem=64GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

# StringTie to assemble transcripts for each sample with the annotation file

F=/data/putnamlab/jillashey/Francois_data/stringTie_hisat2/pdam/BAMfiles/BAMfiles_pdamSamples

array1=($(ls $F/*bam))
for i in ${array1[@]}; do
	stringtie -G /data/putnamlab/jillashey/genome/Pdam/pdam_annotation.gff3 -e -o ${i}.gtf ${i}
	echo "${i}"
done

sbatch stringTie_pdamSamplesOnly_assemble_TuftsGFF.sh
```
Submitted batch job 1672132

b) Merge stringTie gtf results 

```
# Move only pdam samples
mv *gtf ../GTFfiles/GTFfiles_pdamSamples_TuftsGFF

ls *gtf > pdamSamples_TuftsGFF_mergelist.txt
cat pdamSamples_TuftsGFF_mergelist.txt 

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

stringtie --merge -p 8 -G /data/putnamlab/jillashey/genome/Pdam/pdam_annotation.gff3 -o stringtie_pdamSamples_TuftsGFF_merged.gtf pdamSamples_TuftsGFF_mergelist.txt

```

c) Assess assembly quality

```
module load gffcompare/0.11.5-foss-2018b

gffcompare -r /data/putnamlab/jillashey/genome/Pdam/pdam_annotation.gff3 -o ../merged_TuftsGFF stringtie_pdamSamples_TuftsGFF_merged.gtf
```

d) Re-estimate assembly

```
nano stringTie_pdamSamplesOnly_TuftsGFF_re-assemble.sh

#!/bin/bash
#SBATCH --mem=64GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/stringTie_hisat2/pdam/BAMfiles/BAMfiles_pdamSamples

array1=($(ls $F/*bam))
for i in ${array1[@]}
do
stringtie -e -G /data/putnamlab/jillashey/genome/Pdam/pdam_annotation.gff3 -o ${i}.merge.gtf ${i}
echo "${i}"
done

sbatch stringTie_pdamSamplesOnly_TuftsGFF_re-assemble.sh
```

Submitted batch job 1672135

e) Create gene matrix

```
module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/stringTie_hisat2/pdam/GTFfiles_merge/GTFfiles_merge_pdamSamples_TuftsGFF

array2=($(ls *merge.gtf))

for i in ${array2[@]}
do
echo "${i} $F/${i}" >> sample_list_pdamSamples_TuftsGFF.txt
done

python prepDE.py -g gene_count_pdamSamples_TuftsGFF_matrix.csv -i sample_list_pdamSamples_TuftsGFF.txt

```

The csvs were produced, but there are 0 gene counts for every gene in every sample


### Another gff annotation file used below 

StringTie is having problems with the new gff file, so I added gene_id= in the attributes column. Hopefully that will allow it to run correctly. 


a) Assemble reads with genome annotation

```
nano stringTie_pdamSamplesOnly_assemble_rgGFF_fixed.sh

#!/bin/bash
#SBATCH --mem=64GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

# StringTie to assemble transcripts for each sample with the annotation file

F=/data/putnamlab/jillashey/Francois_data/stringTie_hisat2/pdam/BAMfiles/BAMfiles_pdamSamples

array1=($(ls $F/*bam))
for i in ${array1[@]}; do
	stringtie -G /data/putnamlab/jillashey/genome/Pdam/pdam_annotation_Addgene_id_fixed.gtf -e -o ${i}.gtf ${i}
	echo "${i}"
done

sbatch stringTie_pdamSamplesOnly_assemble_rgGFF_fixed.sh
```

Submitted batch job 1672292

b) Merge stringTie gtf results 

```
# Move only pdam samples
mv *gtf ../GTFfiles/GTFfiles_pdamSamples_rgGFF_geneid_fixed

ls *gtf > pdamSamples_rgGFF_geneid_mergelist.txt
cat pdamSamples_rgGFF_geneid_mergelist.txt 

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

stringtie --merge -p 8 -G /data/putnamlab/jillashey/genome/Pdam/pdam_annotation_Addgene_id_fixed.gtf -o stringtie_pdamSamples_rgGFF_geneid_merged.gtf pdamSamples_rgGFF_geneid_mergelist.txt

```

Gave me an output, but also gave me a lot of these: 

```
Warning: invalid start coordinate at line:
###			NA	NA				
Warning: invalid start coordinate at line:
###			NA	NA				
Warning: invalid start coordinate at line:
###			NA	NA				
Warning: invalid start coordinate at line:
###			NA	NA
```

Not sure what that means. Going to continue on, but not sure this gff file is correct

c) Assess assembly quality

```
module load gffcompare/0.11.5-foss-2018b

gffcompare -r /data/putnamlab/jillashey/genome/Pdam/pdam_annotation_Addgene_id_fixed.gtf -o merged_rgGFF_geneid stringtie_pdamSamples_rgGFF_geneid_merged.gtf
```

Gave me the same warnings as above...

I think i was having issues because of the hisat2 index -- it was generated with the NCBI fasta file, not the reef genomics fasta file. Gotta redo hisat2 with reef genomics files. all that will be in other md files 