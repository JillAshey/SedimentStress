### Running HISAT2 and stringtie with pcomp genome

```
cd Francois_data/output/
mkdir hisat2/pcomp/
cd hisat2/pcomp/
ln -s /data/putnamlab/jillashey/Francois_data/data/trimmed/*trim.fq ./
```

Index reference genome 

```
module load HISAT2/2.1.0-foss-2018b 

hisat2-build -f /data/putnamlab/jillashey/genome/Pcomp/Porites_compressa_contigs.fasta ./Pcomp_ref
```

First, try running one sample to genome with hisat2

```
nano Pcomp_hisat2_24_1_test.sh

#!/bin/bash
#SBATCH --mem=10G

module load SAMtools/1.9-foss-2018b
module load HISAT2/2.1.0-foss-2018b 

#Specify working directory
F=/data/putnamlab/jillashey/Francois_data/output/hisat2/plut

# This then makes it into a bam file
# And then also sorts the bam file because Stringtie takes a sorted file for input
# And then removes the sam file because I don't need it anymore

hisat2 -p 8 --dta -q -x Plut_ref -U 24_1.fastq.trim.fq -S 24_1.fastq.trim.fq.sam
samtools sort -@ 8 -o 24_1.fastq.trim.fq.bam 24_1.fastq.trim.fq.sam
   echo "24_1.fastq.trim.fq_bam"
rm 24_1.fastq.trim.fq.sam
   echo "HISAT2 SE 24_1.fastq.trim.fq" $(date)
done

sbatch Pcomp_hisat2_24_1_test.sh

```

Submitted batch job 1669935

Try with all samples now 

```
nano Pcomp_hisat2_allsamples.sh
```

```
#!/bin/bash
#SBATCH --mem=64G

module load SAMtools/1.9-foss-2018b
module load HISAT2/2.1.0-foss-2018b 

#Specify working directory
F=/data/putnamlab/jillashey/Francois_data/output/hisat2/pcomp

array1=($(ls $F/*trim.fq))

# This then makes it into a bam file
# And then also sorts the bam file because Stringtie takes a sorted file for input
# And then removes the sam file because I don't need it anymore

for i in ${array1[@]}; do
        hisat2 -p 8 --dta -q -x Pcomp_ref -U ${i} -S ${i}.sam
        samtools sort -@ 8 -o ${i}.bam ${i}.sam
    		echo "${i}_bam"
        rm ${i}.sam
        echo "HISAT2 SE ${i}" $(date)
done

sbatch Pcomp_hisat2_allsamples.sh

```
Submitted batch job 1671801

