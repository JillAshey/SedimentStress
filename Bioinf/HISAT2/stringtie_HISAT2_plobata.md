StringTie - hisat2, plutea genome

a) Assemble reads with genome annotation

```
nano stringTie_plut_assemble.sh

#!/bin/bash
#SBATCH --mem=64GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

# StringTie to assemble transcripts for each sample with the annotation file

F=/data/putnamlab/jillashey/Francois_data/stringTie_hisat2/plob/BAMfiles/

array1=($(ls $F/*bam))
for i in ${array1[@]}; do
	stringtie -G /data/putnamlab/jillashey/genome/Plutea/Plut.GFFannotation.fixed_transcript.gff -e -o ${i}.gtf ${i}
	echo "${i}"
done

sbatch stringTie_plut_assemble.sh
```
Submitted batch job 1671802

b) Merge stringTie gtf results 

```
mv *gtf ../GTFfiles/

ls *gtf > plob_mergelist.txt
cat plob_mergelist.txt 

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

stringtie --merge -p 8 -G /data/putnamlab/jillashey/genome/Plutea/Plut.GFFannotation.fixed_transcript.gff -o stringtie_plob_merged.gtf plob_mergelist.txt

```

Doing it for only known plob samples: 

```
mkdir GTFfiles_plobSamples
mv 6_1.fastq.trim.fq.bam.gtf 7_1.fastq.trim.fq.bam.gtf 8_1.fastq.trim.fq.bam.gtf 9_1.fastq.trim.fq.bam.gtf 21_1.fastq.trim.fq.bam.gtf 22_1.fastq.trim.fq.bam.gtf 23_1.fastq.trim.fq.bam.gtf 25_1.fastq.trim.fq.bam.gtf 26_1.fastq.trim.fq.bam.gtf 27_1.fastq.trim.fq.bam.gtf 29_1.fastq.trim.fq.bam.gtf 34_1.fastq.trim.fq.bam.gtf GTFfiles_plobSamples/

ls *gtf > plobSamples_mergelist.txt
cat plobSamples_mergelist.txt 

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

stringtie --merge -p 8 -G /data/putnamlab/jillashey/genome/Plutea/Plut.GFFannotation.fixed_transcript.gff -o stringtie_plobSamples_merged.gtf plobSamples_mergelist.txt
```


c) Assess assembly quality

```
module load gffcompare/0.11.5-foss-2018b

gffcompare -r /data/putnamlab/jillashey/genome/Plutea/Plut.GFFannotation.fixed_transcript.gff -o ../merged stringtie_plob_merged.gtf

31126 reference transcripts loaded.
  87046 query transfrags loaded.

```

Doing it for only known plob samples: 

```
module load gffcompare/0.11.5-foss-2018b

gffcompare -r /data/putnamlab/jillashey/genome/Plutea/Plut.GFFannotation.fixed_transcript.gff -o ../plobSamples_merged stringtie_plobSamples_merged.gtf

 31126 reference transcripts loaded.
  53665 query transfrags loaded.
```

To do GOSeq on plob samples only, need merged.annotated file with only those samples. I get that file through the above command


d) Re-estimate assembly

```
nano stringTie_plut_re-assemble.sh

#!/bin/bash
#SBATCH --mem=64GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/stringTie_hisat2/plob/BAMfiles/

array1=($(ls $F/*bam))
for i in ${array1[@]}
do
stringtie -e -G /data/putnamlab/jillashey/genome/Plutea/Plut.GFFannotation.fixed_transcript.gff -o ${i}.merge.gtf ${i}
echo "${i}"
done

sbatch stringTie_plut_re-assemble.sh
```
Submitted batch job 1671843

```
mv *merge.gtf /data/putnamlab/jillashey/Francois_data/stringTie_hisat2/plob/GTFfiles_merged

```

Doing it for only known plob samples: 

```
mv 6_1.fastq.trim.fq.bam 7_1.fastq.trim.fq.bam 8_1.fastq.trim.fq.bam 9_1.fastq.trim.fq.bam 21_1.fastq.trim.fq.bam 22_1.fastq.trim.fq.bam 23_1.fastq.trim.fq.bam 25_1.fastq.trim.fq.bam 26_1.fastq.trim.fq.bam 27_1.fastq.trim.fq.bam 29_1.fastq.trim.fq.bam 34_1.fastq.trim.fq.bam BAMfiles_plobSamples/


nano stringTie_plobSamples_re-assemble.sh

#!/bin/bash
#SBATCH --mem=64GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab

module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/stringTie_hisat2/plob/BAMfiles/BAMfiles_plobSamples/

array1=($(ls $F/*bam))
for i in ${array1[@]}
do
stringtie -e -G /data/putnamlab/jillashey/genome/Plutea/Plut.GFFannotation.fixed_transcript.gff -o ${i}.merge.gtf ${i}
echo "${i}"
done

sbatch stringTie_plobSamples_re-assemble.sh

```
Submitted batch job 1671942


e) Create gene matrix

```
module load StringTie/2.1.1-GCCcore-7.3.0
module load gffcompare/0.11.5-foss-2018b
module load Python/2.7.15-foss-2018b

F=/data/putnamlab/jillashey/Francois_data/stringTie_hisat2/plob/GTFfiles_merged/

array2=($(ls *merge.gtf))

for i in ${array2[@]}
do
echo "${i} $F${i}" >> sample_list_plob.txt
done

python prepDE.py -g gene_count_plob_matrix.csv -i sample_list_plob.txt

```


f) Secure-copy gene counts onto local computer 

```
scp jillashey@bluewaves.uri.edu:/data/putnamlab/jillashey/Francois_data/stringTie_hisat2/plob/GTFfiles_merged/gene_count_plob_matrix.csv /Users/jillashey/Desktop/
```