```
module load BLAST+/2.8.1-foss-2018b 

# protein against protein

# fasta file blasting against - kind of like a reference 
GCF_003704095.1_ASM370409v1_protein.faa

# query fasta file
/data/putnamlab/jillashey/genome/Pdam/ReefGenomics/pdam_proteins.fasta

head ../InterProScan/GCF_003704095.1_ASM370409v1_protein.faa
head /data/putnamlab/jillashey/genome/Pdam/ReefGenomics/pdam_proteins.fasta

# How many sequences are in each file?
grep -c '^>' GCF_003704095.1_ASM370409v1_protein.faa # 25183
grep -c '^>' /data/putnamlab/jillashey/genome/Pdam/ReefGenomics/pdam_proteins.fasta # 26077

# make database to blast against 
makeblastdb -in GCF_003704095.1_ASM370409v1_protein.faa -dbtype prot

# standard
nano RG_vs_NCBI_pdam_protein_blast.sh

#!/bin/bash
#SBATCH --job-name="blastp"
#SBATCH -t 30-00:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="RG_vs_NCBI_pdam_protein_blast_out_error"
#SBATCH --output="RG_vs_NCBI_pdam_protein_blast_out"

module load BLAST+/2.8.1-foss-2018b 

blastp -query /data/putnamlab/jillashey/genome/Pdam/ReefGenomics/pdam_proteins.fasta -db GCF_003704095.1_ASM370409v1_protein.faa -out RG_vs_NCBI_pdam_protein_blast.txt

sbatch RG_vs_NCBI_pdam_protein_blast.sh
```

```
# Statistically sig matches
nano RG_vs_NCBI_pdam_protein_blast.sig.sh

#!/bin/bash
#SBATCH --job-name="blastp"
#SBATCH -t 30-00:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="RG_vs_NCBI_pdam_protein_blast_out_error"
#SBATCH --output="RG_vs_NCBI_pdam_protein_blast_out"

module load BLAST+/2.8.1-foss-2018b 

blastp -query /data/putnamlab/jillashey/genome/Pdam/ReefGenomics/pdam_proteins.fasta \
-db GCF_003704095.1_ASM370409v1_protein.faa \
-out RG_vs_NCBI_pdam_protein_blast.sig.tab \
-evalue 1e-5 \
-outfmt 7

sbatch RG_vs_NCBI_pdam_protein_blast.sig.sh

# Submitted batch job 1713991
```

```
# Best hit only for each query 
nano RG_vs_NCBI_pdam_protein_blast.BestHit.sh

#!/bin/bash
#SBATCH --job-name="blastp"
#SBATCH -t 30-00:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="RG_vs_NCBI_pdam_protein_blast_BestHitout_error"
#SBATCH --output="RG_vs_NCBI_pdam_protein_blast_BestHit_out"

module load BLAST+/2.8.1-foss-2018b 

blastp -query /data/putnamlab/jillashey/genome/Pdam/ReefGenomics/pdam_proteins.fasta \
-db GCF_003704095.1_ASM370409v1_protein.faa \
-out RG_vs_NCBI_pdam_protein_blast.BestHit.tab \
-evalue 1e-5 \
-outfmt 6 \
-max_target_seqs 1

sbatch RG_vs_NCBI_pdam_protein_blast.BestHit.sh

Submitted batch job 1716879
```

```
# uniprot BLAST
# Get uniprot protein list 
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz

makeblastdb -in uniprot_sprot.fasta -dbtype prot -out uniprot_sprot.db

Building a new DB, current time: 09/09/2020 19:41:45
New DB name:   /data/putnamlab/jillashey/annotation/blast/uniprot_sprot.db
New DB title:  uniprot_sprot.fasta
Sequence type: Protein
Keep MBits: T
Maximum file size: 1000000000B
Adding sequences from FASTA; added 563082 sequences in 58.5364 seconds.

# Run BLAST
nano RG_vs_NCBI_pdam_protein_blast.Uniprot.sh

#!/bin/bash
#SBATCH --job-name="blastp"
#SBATCH -t 30-00:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="RG_vs_NCBI_pdam_protein_blast_Uniprot_error"
#SBATCH --output="RG_vs_NCBI_pdam_protein_blast_Uniprot_out"

module load BLAST+/2.8.1-foss-2018b 

blastp -query /data/putnamlab/jillashey/genome/Pdam/ReefGenomics/pdam_proteins.fasta \
-db uniprot_sprot.db \
-out RG_vs_NCBI_pdam_protein_blast.Uniprot.tab \
-evalue 1e-5 \
-outfmt 5 \
-num_threads 96 \
-num_alignments 1

sbatch RG_vs_NCBI_pdam_protein_blast.Uniprot.sh

Submitted batch job 1716880
```

Want to run BLAST on acerv. I want to blast it against a bunch of protein files from different species 

```
module load BLAST+/2.8.1-foss-2018b 

# protein against protein

# query fasta file
/data/putnamlab/jillashey/genome/Acerv/Acerv_assembly_v1.0.protein.fa

# In this run, I'm going to try to blast it remotely. Ie don't make a db and put remote as argument 

# standard
nano acerv_protein_remote_blast.sh

#!/bin/bash
#SBATCH --job-name="blastp"
#SBATCH -t 30-00:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="acerv_protein_remote_blast_out_error"
#SBATCH --output="acerv_protein_remote_blast_out"

module load BLAST+/2.8.1-foss-2018b 

blastp -query /data/putnamlab/jillashey/genome/Acerv/Acerv_assembly_v1.0.protein.fa -remote -out acerv_protein_remote_blast.txt

sbatch acerv_protein_remote_blast.sh
```

Not running, giving me this error: BLAST query/options error: Either a BLAST database or subject sequence(s) must be specified
Please refer to the BLAST+ user manual. so im going to try to 



InterProScan

```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/704/095/GCF_003704095.1_ASM370409v1/GCF_003704095.1_ASM370409v1_protein.faa.gz

nano InterProScan_test.sh
#!/bin/bash
#SBATCH --job-name="InterProScan"
#SBATCH -t 30-00:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu

echo "START $(date)"

# Load module
module load InterProScan/5.44-79.0-foss-2018b 
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

Submitted batch job 1716847
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
