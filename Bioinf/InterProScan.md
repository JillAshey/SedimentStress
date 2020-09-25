InterProScan

Running on the putnam lab node (need to be on uri network to use, am using uri vpn) because it has the newest version of interproscan

```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/704/095/GCF_003704095.1_ASM370409v1/GCF_003704095.1_ASM370409v1_protein.faa.gz

nano InterProScan_test.sh
#!/bin/bash
#SBATCH --job-name="InterProScan"
#SBATCH -t 30-00:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH -p putnamlab
#SBATCH --mem=100GB
#SBATCH --error="acerv_interproscan_out_error"
#SBATCH --output="acerv_interproscan_out"


echo "START $(date)"

# Load module
module load InterProScan/5.46-81.0-foss-2019b
module load Java/11.0.2
java -version

# Run InterProScan
interproscan.sh -version
interproscan.sh -f XML -i GCF_003704095.1_ASM370409v1_protein.faa -b pdam_NCBI.interpro.  -iprlookup -goterms -pa
interproscan.sh -mode convert -f GFF3 -i pdam_NCBI.interpro.xml -b pdam_NCBI.interpro.

# -i is the input data
# -b is the output file base
# -f is formats
# -iprlookup enables mapping
# -goterms is GO Term
# -pa is pathway mapping
# -version displays version number
echo "DONE $(date)"

sbatch InterProScan_test.sh

```

Submitted batch job 10562
Did not work 

Erin code for IPS

```
#!/bin/bash
#SBATCH --job-name="InterProScan"
#SBATCH -t 30-00:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=erin_chille@uri.edu
#SBATCH -p putnamlab

cd /data/putnamlab/erin_chille/mcap2019/annotations/

echo "START $(date)"

# Load module
module load InterProScan/5.46-81.0-foss-2019b
module load Java/11.0.2
java -version

# Run InterProScan
interproscan.sh -version
interproscan.sh -f XML -i ../data/ref/Mcap.IPSprotein.fa -b ./Mcap.interpro.200824  -iprlookup -goterms -pa 
interproscan.sh -mode convert -f GFF3 -i ./Mcap.interpro.200824.xml -b ./Mcap.interpro.200824

# -i is the input data
# -b is the output file base
# -f is formats
# -iprlookup enables mapping
# -goterms is GO Term
# -pa is pathway mapping
# -version displays version number

echo "DONE $(date)"
```

Going to try interproscan on Acerv protein data 

```
nano IPS_acerv.sh

#!/bin/bash
#SBATCH --job-name="InterProScan"
#SBATCH -t 30-00:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu

cd /data/putnamlab/jillashey/annotation/InterProScan/acerv

echo "START $(date)"

# Load module
# module load InterProScan/5.46-81.0-foss-2019b - version erin had in her code, not on bluewaves
module load InterProScan/5.44-79.0-foss-2018b  
module load Java/11.0.2
java -version

# Run InterProScan
interproscan.sh -version
interproscan.sh -f XML -i Acerv_assembly_v1.0.protein.fa -b acerv.interpro -iprlookup -goterms -pa 
interproscan.sh -mode convert -f GFF3 -i acerv.interpro.xml -b acerv.interpro

echo "DONE $(date)"

sbatch IPS_acerv.sh
```

Submitted batch job 1761427 
Error here - there were * in some of the protein sequences and InterProScan aint too happy about that 
Have to remove the * with ```sed```

```
sed -i s/\*//g Acerv_assembly_v1.0.protein.fa

i = in-place (edit file in place)
s = substitute 
/replacement_from_reg_exp/replacement_to_text/ = search and replace statement 
\* = what I want to replace
Add nothing for replacement_to_text
g = global (replace all occurances in file)
```

Submitting acerv IPS job again after removing instances of '*' - Submitted batch job 1761429



Running on the putnam lab node (need to be on uri network to use, am using uri vpn) because it has the newest version of interproscan
