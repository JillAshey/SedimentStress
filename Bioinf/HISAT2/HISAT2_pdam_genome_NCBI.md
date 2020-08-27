### Running HISAT2 and stringtie with pdam genome

```
cd Francois_data/output/
mkdir hisat2/pdam/
cd hisat2/pdam/
ln -s /data/putnamlab/jillashey/Francois_data/data/trimmed/*trim.fq ./
```

Index reference genome 

```
module load HISAT2/2.1.0-foss-2018b 

hisat2-build -f /data/putnamlab/jillashey/genome/Pdam/GCF_003704095.1_ASM370409v1_genomic.fna ./Pdam_ref
```

First, try running one sample to genome with hisat2

```
nano Pdam_hisat2_1_2_test.sh

#!/bin/bash
#SBATCH --mem=10G

module load SAMtools/1.9-foss-2018b
module load HISAT2/2.1.0-foss-2018b 

#Specify working directory
F=/data/putnamlab/jillashey/Francois_data/output/hisat2/pdam

# This then makes it into a bam file
# And then also sorts the bam file because Stringtie takes a sorted file for input
# And then removes the sam file because I don't need it anymore

hisat2 -p 8 --dta -q -x Pdam_ref -U 1_2.fastq.trim.fq -S 1_2.fastq.trim.fq.sam
samtools sort -@ 8 -o 1_2.fastq.trim.fq.bam 1_2.fastq.trim.fq.sam
   echo "1_2.fastq.trim.fq_bam"
rm 1_2.fastq.trim.fq.sam
   echo "HISAT2 SE 1_2.fastq.trim.fq" $(date)
done

sbatch Pdam_hisat2_1_2_test.sh

```
Submitted batch job 1669878
Worked! Got bam file 


Align reads to reference genome 

```
module load SAMtools/1.9-foss-2018b
module load HISAT2/2.1.0-foss-2018b 

nano Pdam_hisat2_allsamples.sh
```

```
#!/bin/bash
#SBATCH --mem=10G

module load SAMtools/1.9-foss-2018b
module load HISAT2/2.1.0-foss-2018b 

#Specify working directory
F=/data/putnamlab/jillashey/Francois_data/output/hisat2/pdam

array1=($(ls $F/*trim.fq))

# Creating sam file, converting to bam file, removing sam file 

for i in ${array1[@]}; do
        hisat2 -p 8 --dta -q -x Pdam_ref -U ${i} -S ${i}.sam
        samtools sort -@ 8 -o ${i}.bam ${i}.sam
    		echo "${i}_bam"
        rm ${i}.sam
        echo "HISAT2 SE ${i}" $(date)
done

```

From erin code, took out --rf

```
chmod u+x Pdam_hisat2_allsamples.sh

# Had some issues with arrays not running all samples together. I think I put too much
# trying with smaller subset first 
sbatch --array 1-5 Pdam_hisat2_allsamples.sh 
Submitted batch job 1669941 
```

Taking quite a bit of time. Slurm outputs say that its running multiple times ? so not sure when it stops. Should ask Erin about this...

This is taking forever and Erin said one sample should only take about 3 min. I'm going to run a single sample again and see how it turns out 


```
nano Pdam_hisat2_28_2_test.sh

#!/bin/bash
#SBATCH --mem=64G

module load SAMtools/1.9-foss-2018b
module load HISAT2/2.1.0-foss-2018b 

#Specify working directory
F=/data/putnamlab/jillashey/Francois_data/output/hisat2/pdam

# This then makes it into a bam file
# And then also sorts the bam file because Stringtie takes a sorted file for input
# And then removes the sam file because I don't need it anymore

hisat2 -p 8 --dta -q -x Pdam_ref -U 28_2.fastq.trim.fq -S 28_2.fastq.trim.fq.sam
samtools sort -@ 8 -o 28_2.fastq.trim.fq.bam 28_2.fastq.trim.fq.sam
   echo "28_2.fastq.trim.fq_bam"
rm 28_2.fastq.trim.fq.sam
   echo "HISAT2 SE 28_2.fastq.trim.fq" $(date)
done

sbatch Pdam_hisat2_28_2_test.sh

```
Submitted batch job 1670724

Based on slurm output, there may be something up with how sam tools is reading the files. Going to try taking out the @ argument in the sam tools command to see if that affects anything and run it on a few

```
sbatch --array 60-64 Pdam_hisat2_allsamples.sh 
Submitted batch job 1670885
```

Let it run its course and finished all of them. Going to run it as just a normal job now to see how long it takes 
Submitted batch job 1671797
 Ran all, most finished but some ran out of memory. 12 samples did not finish so will just run them on their own
 
 1) 4_1
 ```
nano Pdam_hisat2_4_1.sh

#!/bin/bash
#SBATCH --mem=64G

module load SAMtools/1.9-foss-2018b
module load HISAT2/2.1.0-foss-2018b 

#Specify working directory
F=/data/putnamlab/jillashey/Francois_data/output/hisat2/pdam

# This then makes it into a bam file
# And then also sorts the bam file because Stringtie takes a sorted file for input
# And then removes the sam file because I don't need it anymore

hisat2 -p 8 --dta -q -x Pdam_ref -U 4_1.fastq.trim.fq -S 4_1.fastq.trim.fq.sam
samtools sort -o 4_1.fastq.trim.fq.bam 4_1.fastq.trim.fq.sam
   echo "4_1.fastq.trim.fq_bam"
rm 4_1.fastq.trim.fq.sam
   echo "HISAT2 SE 4_1.fastq.trim.fq" $(date)
done

sbatch Pdam_hisat2_4_1.sh

```
Submitted batch job 1671884

 
2) 4_2

 ```
nano Pdam_hisat2_4_2.sh

#!/bin/bash
#SBATCH --mem=64G

module load SAMtools/1.9-foss-2018b
module load HISAT2/2.1.0-foss-2018b 

#Specify working directory
F=/data/putnamlab/jillashey/Francois_data/output/hisat2/pdam

# This then makes it into a bam file
# And then also sorts the bam file because Stringtie takes a sorted file for input
# And then removes the sam file because I don't need it anymore

hisat2 -p 8 --dta -q -x Pdam_ref -U 4_2.fastq.trim.fq -S 4_2.fastq.trim.fq.sam
samtools sort -o 4_2.fastq.trim.fq.bam 4_2.fastq.trim.fq.sam
   echo "4_2.fastq.trim.fq_bam"
rm 4_2.fastq.trim.fq.sam
   echo "HISAT2 SE 4_2.fastq.trim.fq" $(date)
done

sbatch Pdam_hisat2_4_2.sh

```
Submitted batch job 1671885


3) 5_2

 ```
nano Pdam_hisat2_5_2.sh

#!/bin/bash
#SBATCH --mem=64G

module load SAMtools/1.9-foss-2018b
module load HISAT2/2.1.0-foss-2018b 

#Specify working directory
F=/data/putnamlab/jillashey/Francois_data/output/hisat2/pdam

# This then makes it into a bam file
# And then also sorts the bam file because Stringtie takes a sorted file for input
# And then removes the sam file because I don't need it anymore

hisat2 -p 8 --dta -q -x Pdam_ref -U 5_2.fastq.trim.fq -S 5_2.fastq.trim.fq.sam
samtools sort -o 5_2.fastq.trim.fq.bam 5_2.fastq.trim.fq.sam
   echo "5_2.fastq.trim.fq_bam"
rm 4_2.fastq.trim.fq.sam
   echo "HISAT2 SE 5_2.fastq.trim.fq" $(date)
done

sbatch Pdam_hisat2_5_2.sh

```
Submitted batch job 1671886


4) 21_1

 ```
nano Pdam_hisat2_21_1.sh

#!/bin/bash
#SBATCH --mem=64G

module load SAMtools/1.9-foss-2018b
module load HISAT2/2.1.0-foss-2018b 

#Specify working directory
F=/data/putnamlab/jillashey/Francois_data/output/hisat2/pdam

# This then makes it into a bam file
# And then also sorts the bam file because Stringtie takes a sorted file for input
# And then removes the sam file because I don't need it anymore

hisat2 -p 8 --dta -q -x Pdam_ref -U 21_1.fastq.trim.fq -S 21_1.fastq.trim.fq.sam
samtools sort -o 21_1.fastq.trim.fq.bam 21_1.fastq.trim.fq.sam
   echo "21_1.fastq.trim.fq_bam"
rm 21_1.fastq.trim.fq.sam
   echo "HISAT2 SE 21_1.fastq.trim.fq" $(date)
done

sbatch Pdam_hisat2_21_1.sh

```
Submitted batch job 1671887

5) 26_2

 ```
nano Pdam_hisat2_26_2.sh

#!/bin/bash
#SBATCH --mem=64G

module load SAMtools/1.9-foss-2018b
module load HISAT2/2.1.0-foss-2018b 

#Specify working directory
F=/data/putnamlab/jillashey/Francois_data/output/hisat2/pdam

# This then makes it into a bam file
# And then also sorts the bam file because Stringtie takes a sorted file for input
# And then removes the sam file because I don't need it anymore

hisat2 -p 8 --dta -q -x Pdam_ref -U 26_2.fastq.trim.fq -S 26_2.fastq.trim.fq.sam
samtools sort -o 26_2.fastq.trim.fq.bam 26_2.fastq.trim.fq.sam
   echo "26_2.fastq.trim.fq_bam"
rm 26_2.fastq.trim.fq.sam
   echo "HISAT2 SE 26_2.fastq.trim.fq" $(date)
done

sbatch Pdam_hisat2_26_2.sh

```
Submitted batch job 1671888

6) 29_1

 ```
nano Pdam_hisat2_29_1.sh

#!/bin/bash
#SBATCH --mem=64G

module load SAMtools/1.9-foss-2018b
module load HISAT2/2.1.0-foss-2018b 

#Specify working directory
F=/data/putnamlab/jillashey/Francois_data/output/hisat2/pdam

# This then makes it into a bam file
# And then also sorts the bam file because Stringtie takes a sorted file for input
# And then removes the sam file because I don't need it anymore

hisat2 -p 8 --dta -q -x Pdam_ref -U 29_1.fastq.trim.fq -S 29_1.fastq.trim.fq.sam
samtools sort -o 29_1.fastq.trim.fq.bam 29_1.fastq.trim.fq.sam
   echo "29_1.fastq.trim.fq_bam"
rm 29_1.fastq.trim.fq.sam
   echo "HISAT2 SE 29_1.fastq.trim.fq" $(date)
done

sbatch Pdam_hisat2_29_1.sh

```
Submitted batch job 1671889

7) 32_1

 ```
nano Pdam_hisat2_32_1.sh

#!/bin/bash
#SBATCH --mem=64G

module load SAMtools/1.9-foss-2018b
module load HISAT2/2.1.0-foss-2018b 

#Specify working directory
F=/data/putnamlab/jillashey/Francois_data/output/hisat2/pdam

# This then makes it into a bam file
# And then also sorts the bam file because Stringtie takes a sorted file for input
# And then removes the sam file because I don't need it anymore

hisat2 -p 8 --dta -q -x Pdam_ref -U 32_1.fastq.trim.fq -S 32_1.fastq.trim.fq.sam
samtools sort -o 32_1.fastq.trim.fq.bam 32_1.fastq.trim.fq.sam
   echo "32_1.fastq.trim.fq_bam"
rm 32_1.fastq.trim.fq.sam
   echo "HISAT2 SE 32_1.fastq.trim.fq" $(date)
done

sbatch Pdam_hisat2_32_1.sh

```
Submitted batch job 1671890

8) 32_2

 ```
nano Pdam_hisat2_32_2.sh

#!/bin/bash
#SBATCH --mem=64G

module load SAMtools/1.9-foss-2018b
module load HISAT2/2.1.0-foss-2018b 

#Specify working directory
F=/data/putnamlab/jillashey/Francois_data/output/hisat2/pdam

# This then makes it into a bam file
# And then also sorts the bam file because Stringtie takes a sorted file for input
# And then removes the sam file because I don't need it anymore

hisat2 -p 8 --dta -q -x Pdam_ref -U 32_2.fastq.trim.fq -S 32_2.fastq.trim.fq.sam
samtools sort -o 32_2.fastq.trim.fq.bam 32_2.fastq.trim.fq.sam
   echo "32_2.fastq.trim.fq_bam"
rm 32_2.fastq.trim.fq.sam
   echo "HISAT2 SE 32_2.fastq.trim.fq" $(date)
done

sbatch Pdam_hisat2_32_2.sh

```
Submitted batch job 1671891

9) 35_2

 ```
nano Pdam_hisat2_35_2.sh

#!/bin/bash
#SBATCH --mem=64G

module load SAMtools/1.9-foss-2018b
module load HISAT2/2.1.0-foss-2018b 

#Specify working directory
F=/data/putnamlab/jillashey/Francois_data/output/hisat2/pdam

# This then makes it into a bam file
# And then also sorts the bam file because Stringtie takes a sorted file for input
# And then removes the sam file because I don't need it anymore

hisat2 -p 8 --dta -q -x Pdam_ref -U 35_2.fastq.trim.fq -S 35_2.fastq.trim.fq.sam
samtools sort -o 35_2.fastq.trim.fq.bam 35_2.fastq.trim.fq.sam
   echo "35_2.fastq.trim.fq_bam"
rm 35_2.fastq.trim.fq.sam
   echo "HISAT2 SE 35_2.fastq.trim.fq" $(date)
done

sbatch Pdam_hisat2_35_2.sh

```
Submitted batch job 1671892


10) 36_2

 ```
nano Pdam_hisat2_36_2.sh

#!/bin/bash
#SBATCH --mem=64G

module load SAMtools/1.9-foss-2018b
module load HISAT2/2.1.0-foss-2018b 

#Specify working directory
F=/data/putnamlab/jillashey/Francois_data/output/hisat2/pdam

# This then makes it into a bam file
# And then also sorts the bam file because Stringtie takes a sorted file for input
# And then removes the sam file because I don't need it anymore

hisat2 -p 8 --dta -q -x Pdam_ref -U 36_2.fastq.trim.fq -S 36_2.fastq.trim.fq.sam
samtools sort -o 36_2.fastq.trim.fq.bam 36_2.fastq.trim.fq.sam
   echo "36_2.fastq.trim.fq_bam"
rm 36_2.fastq.trim.fq.sam
   echo "HISAT2 SE 36_2.fastq.trim.fq" $(date)
done

sbatch Pdam_hisat2_36_2.sh

```
Submitted batch job 1671893


11) 37_2

 ```
nano Pdam_hisat2_37_2.sh

#!/bin/bash
#SBATCH --mem=64G

module load SAMtools/1.9-foss-2018b
module load HISAT2/2.1.0-foss-2018b 

#Specify working directory
F=/data/putnamlab/jillashey/Francois_data/output/hisat2/pdam

# This then makes it into a bam file
# And then also sorts the bam file because Stringtie takes a sorted file for input
# And then removes the sam file because I don't need it anymore

hisat2 -p 8 --dta -q -x Pdam_ref -U 37_2.fastq.trim.fq -S 37_2.fastq.trim.fq.sam
samtools sort -o 37_2.fastq.trim.fq.bam 37_2.fastq.trim.fq.sam
   echo "37_2.fastq.trim.fq_bam"
rm 37_2.fastq.trim.fq.sam
   echo "HISAT2 SE 37_2.fastq.trim.fq" $(date)
done

sbatch Pdam_hisat2_37_2.sh

```
Submitted batch job 1671894

12) 39_1

 ```
nano Pdam_hisat2_39_1.sh

#!/bin/bash
#SBATCH --mem=64G

module load SAMtools/1.9-foss-2018b
module load HISAT2/2.1.0-foss-2018b 

#Specify working directory
F=/data/putnamlab/jillashey/Francois_data/output/hisat2/pdam

# This then makes it into a bam file
# And then also sorts the bam file because Stringtie takes a sorted file for input
# And then removes the sam file because I don't need it anymore

hisat2 -p 8 --dta -q -x Pdam_ref -U 39_1.fastq.trim.fq -S 39_1.fastq.trim.fq.sam
samtools sort -o 39_1.fastq.trim.fq.bam 39_1.fastq.trim.fq.sam
   echo "39_1.fastq.trim.fq_bam"
rm 39_1.fastq.trim.fq.sam
   echo "HISAT2 SE 39_1.fastq.trim.fq" $(date)
done

sbatch Pdam_hisat2_39_1.sh
```

Submitted batch job 1671895

### Running HISAT2 on pdam with ref annotation file and fasta file from Reef Genomics

Tufts people were using different genomic resources for pdam (from Reef Genomics), whereas I got my genomic resources from NCBI. The resources from Reef Genomics are much better annotated (includes GO terms). So, I'm going to do hisat2 with fasta file from reef genomics 

```
cd Francois_data/output/
mkdir hisat2/pdam/
cd hisat2/pdam/
ln -s /data/putnamlab/jillashey/Francois_data/data/trimmed/*trim.fq ./
```

Index reference genome 

```
module load HISAT2/2.1.0-foss-2018b 

hisat2-build -f /data/putnamlab/jillashey/genome/Pdam/pdam_scaffolds.fasta ./Pdam_rg_ref
```

First, try running one sample to pdam reef genomics genome with hisat2

```
nano Pdam_hisat2_11_2_test.sh

#!/bin/bash
#SBATCH --mem=64G

module load SAMtools/1.9-foss-2018b
module load HISAT2/2.1.0-foss-2018b 

#Specify working directory
F=/data/putnamlab/jillashey/Francois_data/output/hisat2/pdam_rg

# This then makes it into a bam file
# And then also sorts the bam file because Stringtie takes a sorted file for input
# And then removes the sam file because I don't need it anymore

hisat2 -p 8 --dta -q -x Pdam_rg_ref -U 11_2.fastq.trim.fq -S 11_2.fastq.trim.fq.sam
samtools sort -o 11_2.fastq.trim.fq.bam 11_2.fastq.trim.fq.sam
   echo "11_2.fastq.trim.fq_bam"
rm 11_2.fastq.trim.fq.sam
   echo "HISAT2 SE 11_2.fastq.trim.fq" $(date)
done

sbatch Pdam_hisat2_11_2_test.sh
```
Submitted batch job 1672530
Cool that worked. Now time to run on all 64 samples 

```
nano Pdam_hisat2_allsamples_rgGFF.sh

#!/bin/bash
#SBATCH --mem=100G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu

module load SAMtools/1.9-foss-2018b
module load HISAT2/2.1.0-foss-2018b 

#Specify working directory
F=/data/putnamlab/jillashey/Francois_data/output/hisat2/pdam_rg

array1=($(ls $F/*trim.fq))

# Creating sam file, converting to bam file, removing sam file 

for i in ${array1[@]}; do
        hisat2 -p 8 --dta -q -x Pdam_rg_ref -U ${i} -S ${i}.sam
        samtools sort -@ 8 -o ${i}.bam ${i}.sam
    		echo "${i}_bam"
        rm ${i}.sam
        echo "HISAT2 SE ${i}" $(date)
done

sbatch Pdam_hisat2_allsamples_rgGFF.sh
```
Submitted batch job 1672534





