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



Running on the putnam lab node (need to be on uri network to use, am using uri vpn) because it has the newest version of interproscan. Used IPS_acerv.sh code above. 

Going to try running InterProScan on Ofav now. IPS seems to take a while, like weeks. 

```
nano IPS_ofav.sh

#!/bin/bash
#SBATCH --job-name="IPS_ofav"
#SBATCH -t 30-00:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH -p putnamlab
#SBATCH --mem=200GB
#SBATCH --error="ofav_interproscan_out_error"
#SBATCH --output="ofav_interproscan_out"

# must run on Putnam lab node

cd /data/putnamlab/jillashey/annotation/InterProScan/ofav

echo "START $(date)"

# Load module
module load InterProScan/5.46-81.0-foss-2019b
module load Java/11.0.2
java -version

# Run InterProScan
interproscan.sh -version
interproscan.sh -f XML -i GCF_002042975.1_ofav_dov_v1_protein.faa -b ofav.interpro -iprlookup -goterms -pa 
interproscan.sh -mode convert -f GFF3 -i ofav.interpro.xml -b ofav.interpro

echo "DONE $(date)"

sbatch IPS_ofav.sh
```
Submitted batch job 10661

Ofav finished super fast! I think it was because ofav info is on NCBI so quicker

Going to try running InterProScan on Mcav now. IPS seems to take a while, like weeks. but it did finish super fast on Ofav 

```
nano IPS_mcav.sh

#!/bin/bash
#SBATCH --job-name="IPS_mcav"
#SBATCH -t 30-00:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH -p putnamlab
#SBATCH --mem=200GB
#SBATCH --error="mcav_interproscan_out_error"
#SBATCH --output="mcav_interproscan_out"

# must run on Putnam lab node

cd /data/putnamlab/jillashey/annotation/InterProScan/mcav

echo "START $(date)"

# Load module
module load InterProScan/5.46-81.0-foss-2019b
module load Java/11.0.2
java -version

# Run InterProScan
interproscan.sh -version
interproscan.sh -f XML -i Mcavernosa.maker.proteins.fasta -b mcav.interpro -iprlookup -goterms -pa 
interproscan.sh -mode convert -f GFF3 -i mcav.interpro.xml -b mcav.interpro

echo "DONE $(date)"

sbatch IPS_mcav.sh
```
Submitted batch job 10907

Attempted IPS with pdam, didn't work, going to try again

```
nano IPS_pdam.sh

#!/bin/bash
#SBATCH --job-name="IPS_pdam"
#SBATCH -t 30-00:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH -p putnamlab
#SBATCH --mem=200GB
#SBATCH --error="pdam_interproscan_out_error"
#SBATCH --output="pdam_interproscan_out"

# must run on Putnam lab node

cd /data/putnamlab/jillashey/annotation/InterProScan/pdam

echo "START $(date)"

# Load module
module load InterProScan/5.46-81.0-foss-2019b
module load Java/11.0.2
java -version

# Run InterProScan
interproscan.sh -version
interproscan.sh -f XML -i GCF_003704095.1_ASM370409v1_protein.faa -b pdam.interpro -iprlookup -goterms -pa 
interproscan.sh -mode convert -f GFF3 -i pdam.interpro.xml -b pdam.interpro

echo "DONE $(date)"

sbatch IPS_pdam.sh
```
Submitted batch job 11232

Trying IPS with plob now

```
nano IPS_plob.sh

#!/bin/bash
#SBATCH --job-name="IPS_plob"
#SBATCH -t 30-00:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH -p putnamlab
#SBATCH --mem=200GB
#SBATCH --error="plob_interproscan_out_error"
#SBATCH --output="plob_interproscan_out"

# must run on Putnam lab node

cd /data/putnamlab/jillashey/annotation/InterProScan/plob

echo "START $(date)"

# Load module
module load InterProScan/5.46-81.0-foss-2019b
module load Java/11.0.2
java -version

# Run InterProScan
interproscan.sh -version
interproscan.sh -f XML -i plut2v1.1.proteins.fasta -b plob.interpro -iprlookup -goterms -pa 
interproscan.sh -mode convert -f GFF3 -i plob.interpro.xml -b plob.interpro

echo "DONE $(date)"

sbatch IPS_plob.sh
```

Submitted batch job 11235
