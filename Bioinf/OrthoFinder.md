## Orthofinder 

[OrthoFinder github](https://github.com/davidemms/OrthoFinder)

[OrthoFinder manual](https://github.com/davidemms/OrthoFinder/blob/master/OrthoFinder-manual.pdf)

[OrthoFinder tutorials](https://davidemms.github.io/menu/about.html)

###Running OrthoFinder 

*I did this anaylsis on Bluewaves*

Make new directories 

```
cd /data/putnamlab/jillashey
mkdir mkdir OrthoFinder
cd mkdir OrthoFinder
```

Put protein sequences of interest for all species in folder together 

```
mkdir protein_seqs
cd protein_seqs

# Make softlink so don't have to make copies of protein seqs and move them into new folder 
ln -s /data/putnamlab/jillashey/genome/Acerv/Acerv_assembly_v1.0.protein.fa
ln -s /data/putnamlab/jillashey/genome/Mcav/Mcavernosa_annotation/Mcavernosa.maker.proteins.fasta
ln -s /data/putnamlab/jillashey/genome/Ofav/GCF_002042975.1_ofav_dov_v1_protein.faa
ln -s /data/putnamlab/jillashey/genome/Pacuta/Pocillopora_acuta_HIv1.genes.pep.faa 
ln -s /data/putnamlab/jillashey/genome/Plutea/plut2v1.1.proteins.fasta
```

Run OrthoFinder 

Code from HP Mcap/Pacuta

```
nano orthofinder.sh 

#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=100GB
#SBATCH --export=NONE
#SBATCH --account=putnamlab
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -D /data/putnamlab/jillashey/OrthoFinder
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --error="OrthoFinder_out_error"
#SBATCH --output="OrthoFinder_out"

# load modules needed
# make sure that the foss-year all match

module load DIAMOND/0.9.22-foss-2018b  
module load MCL/14.137-foss-2018b
module load FastME/2.1.6.1-foss-2018b
module load BLAST+/2.8.1-foss-2018b
module load OrthoFinder/2.3.3-foss-2018b-Python-2.7.15

# orthofinder runs as a python script so have to specify .py at the end of it
# requires path to Fastas directory
# using 10 threads, matching the SLRUM parameters above

orthofinder.py -f protein_seqs/ -t 10

sbatch orthofinder.sh 
```

Submitted batch job 1963934

Lets see if this works. Worked! 

Lots of output files/directories

```
cd /data/putnamlab/jillashey/OrthoFinder/protein_seqs/OrthoFinder/Results_Jan25

```

Under ```/data/putnamlab/jillashey/OrthoFinder/protein_seqs/OrthoFinder/Results_Jan25/Comparative_Genomics_Statistics/Statistics_Overall.tsv```:


```
Number of species       5
Number of genes 161090
Number of genes in orthogroups  131269
Number of unassigned genes      29821
Percentage of genes in orthogroups      81.5
Percentage of unassigned genes  18.5
Number of orthogroups   19434
Number of species-specific orthogroups  213
Number of genes in species-specific orthogroups 1010
Percentage of genes in species-specific orthogroups     0.6
Mean orthogroup size    6.8
Median orthogroup size  5.0
G50 (assigned genes)    7
G50 (all genes) 6
O50 (assigned genes)    5124
O50 (all genes) 7441
Number of orthogroups with all species present  9658
Number of single-copy orthogroups       3173
Date    2022-01-25
Orthogroups file        Orthogroups.tsv
Unassigned genes file   Orthogroups_UnassignedGenes.tsv
Per-species statistics  Statistics_PerSpecies.tsv
Overall statistics      Statistics_Overall.tsv
Orthogroups shared between species      Orthogroups_SpeciesOverlaps.tsv

Average number of genes per-species in orthogroup       Number of orthogroups   Percentage of orthogroups       Number of genes Percentage of genes
<1      5645    29.0    17207   13.1
'1      11215   57.7    68787   52.4
'2      1502    7.7     17254   13.1
'3      464     2.4     7678    5.8
'4      215     1.1     4634    3.5
'5      140     0.7     3743    2.9
'6      65      0.3     2074    1.6
'7      60      0.3     2199    1.7
'8      32      0.2     1352    1.0
'9      27      0.1     1266    1.0
'10     10      0.1     528     0.4
11-15   43      0.2     2782    2.1
16-20   8       0.0     722     0.6
21-50   8       0.0     1043    0.8
51-100  0       0.0     0       0.0
101-150 0       0.0     0       0.0
151-200 0       0.0     0       0.0
201-500 0       0.0     0       0.0
501-1000        0       0.0     0       0.0
'1001+  0       0.0     0       0.0

Number of species in orthogroup Number of orthogroups
1       213
2       3155
3       2489
4       3919
5       9658
```

So much info...I need to find out what files will be of use 

Secure-copy files onto computer 

```
cd /Users/jillashey/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/OrthoFinder

scp jillashey@bluewaves.uri.edu:/data/putnamlab/jillashey/OrthoFinder/protein_seqs/OrthoFinder/Results_Jan25/Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv .

scp jillashey@bluewaves.uri.edu:/data/putnamlab/jillashey/OrthoFinder/protein_seqs/OrthoFinder/Results_Jan25/Species_Tree/SpeciesTree_rooted.txt .

scp -r jillashey@bluewaves.uri.edu:/data/putnamlab/jillashey/OrthoFinder/protein_seqs/OrthoFinder/Results_Jan25/Orthologues/Orthologues_Acerv_assembly_v1.0.protein/ .

scp -r jillashey@bluewaves.uri.edu:/data/putnamlab/jillashey/OrthoFinder/protein_seqs/OrthoFinder/Results_Jan25/Orthologues/Orthologues_GCF_002042975.1_ofav_dov_v1_protein .

scp -r jillashey@bluewaves.uri.edu:/data/putnamlab/jillashey/OrthoFinder/protein_seqs/OrthoFinder/Results_Jan25/Orthologues/Orthologues_Mcavernosa.maker.proteins .

scp -r jillashey@bluewaves.uri.edu:/data/putnamlab/jillashey/OrthoFinder/protein_seqs/OrthoFinder/Results_Jan25/Orthologues/Orthologues_plut2v1.1.proteins .

scp -r jillashey@bluewaves.uri.edu:/data/putnamlab/jillashey/OrthoFinder/protein_seqs/OrthoFinder/Results_Jan25/Orthologues/Orthologues_Pocillopora_acuta_HIv1.genes.pep .

scp jillashey@bluewaves.uri.edu:/data/putnamlab/jillashey/OrthoFinder/protein_seqs/OrthoFinder/Results_Jan25/Gene_Trees/OG0000000_tree.txt .

scp jillashey@bluewaves.uri.edu:/data/putnamlab/jillashey/OrthoFinder/protein_seqs/OrthoFinder/Results_Jan25/Gene_Duplication_Events/SpeciesTree_Gene_Duplications_0.5_Support.txt .


```