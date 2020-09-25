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
Please refer to the BLAST+ user manual. 

BLAST with references. Going to use a few different protein refs from different species to see if I can compare acerv proteins to proteins in all other coral species at once 

```
# make database to blast against 
makeblastdb -in /data/putnamlab/jillashey/genome/Adig/GCF_000222465.1_Adig_1.1_protein.faa /data/putnamlab/jillashey/genome/Ofav/GCF_002042975.1_ofav_dov_v1_protein.faa -dbtype prot

Nope cant do that. can only do one at a time. Taking out ofav and moving adig protein file here
ln -s /data/putnamlab/jillashey/genome/Adig/GCF_000222465.1_Adig_1.1_protein.faa .

makeblastdb -in GCF_000222465.1_Adig_1.1_protein.faa -dbtype prot

nano acerv_adig_blastp.sh

#!/bin/bash
#SBATCH --job-name="blastp"
#SBATCH -t 30-00:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="acerv_adig_blastp_out_error"
#SBATCH --output="acerv_adig_blastp_out"

module load BLAST+/2.8.1-foss-2018b 

# Comparing Acerv proteins against Adig proteins 

blastp -query Acerv_assembly_v1.0.protein.fa -db GCF_000222465.1_Adig_1.1_protein.faa -out acerv_adig_blastp.sig.txt -evalue 1e-5 -outfmt 7

sbatch acerv_adig_blastp.sh

```
Submitted batch job 1761749

Worked! Maybe try running against a different species? Let's try Ofav

```
ln -s /data/putnamlab/jillashey/genome/Ofav/GCF_002042975.1_ofav_dov_v1_protein.faa .

makeblastdb -in GCF_002042975.1_ofav_dov_v1_protein.faa -dbtype prot

Illegal instruction (core dumped) -- gives me this when I try to make the ofav protein db. Weird. Guess I need to try another species 
```

BLAST is only helpful if the protein you are blasting against are labelled already with protein names

Let's try with



Want to try running Diamond Blast 

```
nano acerv_diamond_blastp.sh

#!/bin/bash
#SBATCH --job-name="diamond-blastp"
#SBATCH -t 30-00:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --mem=100GB
#SBATCH --error="acerv_diamond_blastp_out_error"
#SBATCH --output="acerv_diamond_blastp_out"

echo "START" $(date)
module load DIAMOND/2.0.0-GCC-8.3.0 #Load DIAMOND

diamond blastp -d /data/putnamlab/shared/databases/nr.dmnd -q Acerv_assembly_v1.0.protein.fa -o acerv_annot -f 100 -b20 --more-sensitive -e 0.00001 -k1

echo "Search complete... converting format to XML and tab"

diamond view -a acerv_annot.daa -o acerv_annot.xml -f 5
diamond view -a acerv_annot.daa -o acerv_annot.tab -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen

echo "STOP" $(date)

sbatch acerv_diamond_blastp.sh

```
Submitted batch job 1761752

Maybe I could do more if I binded protein files of different species together???
