StringTie for pdam genome - Francois 64 samples - with STAR

Copy BAM files to new stringTie directory 

a) Assemble reads with genome annotation

```
nano stringTie_assemble_pdam.sh
```

```
#!/bin/bash
#SBATCH -t 48:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Francois_data/scripts
#SBATCH --error="stringTie_out_error"
#SBATCH --output="stringTie_out"

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

# StringTie to assemble transcripts for each sample with the annotation file

F=/data/putnamlab/jillashey/Francois_data/stringTie/pdam/BAMfiles

array1=($(ls $F/*Aligned.sortedByCoord.out.bam))
for i in ${array1[@]}
do
stringtie -G /data/putnamlab/jillashey/genome/Pdam/GCF_003704095.1_ASM370409v1_genomic.gff -o ${i}.gtf ${i}
echo "${i}"
done
```

```
sbatch stringTie_assemble_pdam.sh
```
Submitted batch job 1667439


b) Merge stringTie GTF files 

Attempt 1

```
mv *Aligned.sortedByCoord.out.bam.gtf /data/putnamlab/jillashey/Francois_data/stringTie/pdam/GTFfiles/

ls *gtf > pdam_mergelist.txt
cat pdam_mergelist.txt

sed 's,^,'$PWD' ,' pdam_mergelist.txt # add path in front of gtf files

# There is a space in between sampleID and working directory
# I am a cheater and did it in R. Secure copied pdam_mergelist.txt to local computer, removed space between path and sample name in R, secure copied results to server as pdam_mergelist_no_spaces.txt
pdam_mergelist_no_spaces.txt

sed '1d' pdam_mergelist_no_spaces.txt > test_remove_firstrow.txt # remove header row that was put there in R
```

```
# Dont need to submit as job

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

stringtie --merge -p 8 -G /data/putnamlab/jillashey/genome/Pdam/GCF_003704095.1_ASM370409v1_genomic.gff -o stringtie_pdam_merged.gtf /data/putnamlab/jillashey/Francois_data/stringTie/pdam/GTFfiles/test_remove_firstrow.txt
```

c) Assess quality of assembly

```
module load gffcompare/0.11.5-foss-2018b

gffcompare -r /data/putnamlab/jillashey/genome/Pdam/GCF_003704095.1_ASM370409v1_genomic.gff -o ../merged /data/putnamlab/jillashey/Francois_data/stringTie/pdam/GTFfiles/stringtie_pdam_merged.gtf
```

d) Re-estimate assembly / abundance 

```
nano stringTie_reassemble_pdam.sh
```

```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=50GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Francois_data/scripts
#SBATCH --error="stringTie_reassemble_out_error"
#SBATCH --output="stringTie_reassemble_out"

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/stringTie/pdam/BAMfiles

array1=($(ls $F/*Aligned.sortedByCoord.out.bam))

for i in ${array1[@]}
do
stringtie -e -G /data/putnamlab/jillashey/Francois_data/stringTie/pdam/GTFfiles/stringtie_pdam_merged.gtf -o ${i}.merge.gtf ${i}
echo "${i}"
done 	 
```

```
sbatch stringTie_reassemble_pdam.sh
```

Submitted batch job 1667634

e) Create gene matrix

```
module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/stringTie/pdam/GTFfiles_merged/

array2=($(ls *merge.gtf))

for i in ${array2[@]}
do
	echo "${i} $F/${i}" >> sample_list.txt
done

awk '{print $2}' sample_list.txt >> sample_list_pdam.txt

#for i in ${array2[@]}
#do
#echo "$(echo ${i}|sed 's,^,'$PWD' ,') $F/${i}" >> sample_list_pdam.txt
#done


# do python step in directory where prepDE.py code is
python prepDE.py -g gene_count_pdam_matrix.csv -i /data/putnamlab/jillashey/Francois_data/stringTie/pdam/GTFfiles_merged/sample_list_full_paths_pdam.txt

python prepDE.py -g gene_count_pdam_matrix.csv -i sample_list_pdam.txt
```

Not working

Error: line should have a sample ID and a file path:
/data/putnamlab/jillashey/Francois_data/stringTie/pdam/GTFfiles_merged/10_2.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf


```
F=/data/putnamlab/jillashey/Francois_data/stringTie/pdam/GTFfiles_merged/

array2=($(ls *merge.gtf))

for i in ${array2[@]}
do
echo -e “${i} $F/${i} \r\n” >> sample_list.txt
done

python prepDE.py -g gene_count_pdam_matrix.csv -i /data/putnamlab/jillashey/Francois_data/stringTie/pdam/GTFfiles_merged/sample_list.txt


path=/data/putnamlab/jillashey/Francois_data/stringTie/pdam/GTFfiles_merged/

for $sample in “$path”/*.”$.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf”
do

test_path = {$file //“$path””10_2.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf”}

sample_id= {10_2.fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf//$end} ```

Straight from my connelly script that worked

```
module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/stringTie/pdam/GTFfiles_merged/

array2=($(ls *merge.gtf))

for i in ${array2[@]}
do
	echo "${i} $F${i}" >> sample_list_connelly_stringTie_example_script.txt
done

python prepDE.py -g gene_count_pdam_matrix.csv -i sample_list_connelly_stringTie_example_script.txt


array2=($(ls *merge.gtf))

for i in ${array2[@]}
do
	echo "${i} $F" >> sample_list_connelly_stringTie_example_script_pleasework.txt
done

python prepDE.py -g gene_count_pdam_matrix.csv -i sample_list_connelly_stringTie_example_script_pleasework.txt

python prepDE.py -g gene_count_pdam_matrix.csv -i sample_list_connelly_stringTie_example_script.txt
```


### Fixed GFF file with GO terms

I recently did some fancy manuvering and adding GO terms to the NCBI gff file (it had none before, only the reef genomics file included them, but the RG file had some discrepancies in gene counts). The R script for the GO-term additions is [here](https://github.com/JillAshey/SedimentStress/blob/master/RAnalysis/AddAnnotations_NCBI_pdam.R). Now I'm going to use the fixed gff file to run STAR with pdam. 

a) Assemble reads with genome annotation

```
nano stringTie_assemble_pdam_GOterms.sh

#!/bin/bash
#SBATCH -t 48:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="stringTie_pdam_GOterms_out_error"
#SBATCH --output="stringTie_pdam_GOterms_out"

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

# StringTie to assemble transcripts for each sample with the annotation file

F=/data/putnamlab/jillashey/Francois_data/Hawaii/stringTie_star/pdam_GOterms/BAM

array1=($(ls $F/*Aligned.sortedByCoord.out.bam))
for i in ${array1[@]}
do
stringtie -G /data/putnamlab/jillashey/genome/Pdam/NCBI/pdam_NCBI_annotation_fixed_GOterms.gff -o ${i}.gtf ${i}
echo "${i}"
done

sbatch stringTie_assemble_pdam_GOterms.sh
``` 

Submitted batch job 1701979

b) Merge stringTie gtf results 

```
mv *gtf ../GTFfiles/

ls *gtf > pdam_GOterms_mergelist.txt
cat pdam_GOterms_mergelist.txt

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

stringtie --merge -p 8 -G /data/putnamlab/jillashey/genome/Pdam/NCBI/pdam_NCBI_annotation_fixed_GOterms.gff -o stringtie_pdam_GOterms_merged.gtf pdam_GOterms_mergelist.txt
```

Lots of errors that look like this: 

```
Error: discarding overlapping duplicate gene feature (378199-378285) with ID=gene-Trnay-gua-10
Error: discarding overlapping duplicate gene feature (378199-378285) with ID=gene-Trnay-gua-10
Error: discarding overlapping duplicate gene feature (378199-378285) with ID=gene-Trnay-gua-10
Error: discarding overlapping duplicate gene feature (5307-5379) with ID=gene-Trnay-gua-4
Error: discarding overlapping duplicate gene feature (5307-5379) with ID=gene-Trnay-gua-4
Error: discarding overlapping duplicate gene feature (5307-5379) with ID=gene-Trnay-gua-4
Error: discarding overlapping duplicate gene feature (5307-5379) with ID=gene-Trnay-gua-4
Error: discarding overlapping duplicate gene feature (5307-5379) with ID=gene-Trnay-gua-4
Error: discarding overlapping duplicate gene feature (5307-5379) with ID=gene-Trnay-gua-4
Error: discarding overlapping duplicate gene feature (5307-5379) with ID=gene-Trnay-gua-4
Error: discarding overlapping duplicate gene feature (5307-5379) with ID=gene-Trnay-gua-4
Error: discarding overlapping duplicate gene feature (5307-5379) with ID=gene-Trnay-gua-4
Error: discarding overlapping duplicate gene feature (5307-5379) with ID=gene-Trnay-gua-4
Error: discarding overlapping duplicate gene feature (5307-5379) with ID=gene-Trnay-gua-4
```

c) Assess assembly quality

```
module load gffcompare/0.11.5-foss-2018b

gffcompare -r /data/putnamlab/jillashey/genome/Pdam/NCBI/pdam_NCBI_annotation_fixed_GOterms.gff -o pdam_GOterms.merged stringtie_pdam_GOterms_merged.gtf

```

Gave me these errors: 

```
Error: discarding overlapping duplicate transcript feature (269346-269419) with ID=rna-Trnat-agu-15
Error: discarding overlapping duplicate transcript feature (269346-269419) with ID=rna-Trnat-agu-15
Error: discarding overlapping duplicate transcript feature (269346-269419) with ID=rna-Trnat-agu-15
Error: discarding overlapping duplicate transcript feature (269346-269419) with ID=rna-Trnat-agu-15
Error: discarding overlapping duplicate transcript feature (269346-269419) with ID=rna-Trnat-agu-15
Error: discarding overlapping duplicate transcript feature (269346-269419) with ID=rna-Trnat-agu-15
Error: discarding overlapping duplicate transcript feature (269346-269419) with ID=rna-Trnat-agu-15
Error: discarding overlapping duplicate transcript feature (269346-269419) with ID=rna-Trnat-agu-15
Error: discarding overlapping duplicate transcript feature (269346-269419) with ID=rna-Trnat-agu-15
Error: discarding overlapping duplicate transcript feature (269346-269419) with ID=rna-Trnat-agu-15
```

d) Re-estimate assembly 

```
nano stringTie_assemble_pdam_GOterms_re-assemble.sh

#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="re-assemble_pdam_GOterms_out_error"
#SBATCH --output="re-assemble_pdam_GOterms_out"

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/Hawaii/stringTie_star/pdam_GOterms/BAM

array1=($(ls $F/*bam))
for i in ${array1[@]}
do
stringtie -e -G /data/putnamlab/jillashey/genome/Pdam/NCBI/pdam_NCBI_annotation_fixed_GOterms.gff -o ${i}.merge.gtf ${i}
echo "${i}"
done

sbatch stringTie_assemble_pdam_GOterms_re-assemble.sh

mv *merge.gtf ../GTF_merge
```

Submitted batch job 1702705

e) Create gene matrix

```
module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/Hawaii/stringTie_star/pdam_GOterms/GTF_merge/

array2=($(ls *merge.gtf))

for i in ${array2[@]}
do
echo "${i} $F${i}" >> sample_list_pdam_GOterms.txt
done

python prepDE.py -g gene_count_pdam_GOterms_matrix.csv -i sample_list_pdam_GOterms.txt
```

f) Secure-copy gene counts onto local computer

```
scp jillashey@bluewaves.uri.edu:/data/putnamlab/jillashey/Francois_data/Hawaii/stringTie_star/pdam_GOterms/GTF_merge/gene_count_pdam_GOterms_matrix.csv /Users/jillashey/Desktop/Putnamlab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/star/
```