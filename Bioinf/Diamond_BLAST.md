## Diamond BLAST 

Go to the sbatch_executables subdirectory in the Putnam Lab shared folder and run the script, make_diamond_nr_db.sh. This script, created by Erin Chille on August 6, 2020, downloads the most recent nr database in FASTA format from NCBI and uses it to make a Diamond-formatted nr database.

```
#!/bin/bash
#SBATCH --job-name="diamond_db_update" #CHANGE_NAME
#SBATCH -t 100:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu #CHANGE_EMAIL
#SBATCH --mem=128GB

echo "START" $(date)
module load DIAMOND/2.0.0-GCC-8.3.0 #Load DIAMOND

cd /data/putnamlab/shared/databases/

echo "Downloading nr database" $(date)
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz #download nr database in fasta format

echo "Making diamond database" $date
diamond makedb --in nr.gz -d nr
diamond dbinfo -d nr.dmnd

echo "STOP" $(date)
```

Submitted batch job 13094

Previously, I ran blastp with protein queries against protein db, but here I'm going to try blastx with transcript queries against protein db. 

*Run on Putnam Lab Node*

### Florida species

#### Acerv 

Ran previously, failed before I got any results due to exceeded memory. Rerun on Putnam lab node

```
# Check number of genes 
zgrep -c "^>" Acerv_assembly_v1.0.mRNA.fa
33322 

nano acerv_diamond_blastx.sh

#!/bin/bash
#SBATCH --job-name="diamond-blastx"
#SBATCH -t 240:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH -p putnamlab
#SBATCH --mem=100GB
#SBATCH --error="acerv_diamond_blastx_out_error"
#SBATCH --output="acerv_diamond_blastx_out"

echo "START" $(date)
module load DIAMOND/2.0.0-GCC-8.3.0 #Load DIAMOND

echo "Updating Acerv annotation" $(date)
diamond blastx -d /data/putnamlab/shared/databases/nr.dmnd -q Acerv_assembly_v1.0.mRNA.fa -o Acerv_annot -f 100 -b20 --more-sensitive -e 0.00001 -k1

echo "Search complete... converting format to XML and tab"

diamond view -a Acerv_annot.daa -o Acerv_annot.xml -f 5
diamond view -a Acerv_annot.daa -o Acerv_annot.tab -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen

echo "STOP" $(date)

sbatch acerv_diamond_blastx.sh

```
Submitted batch job 13114


#### Mcav 

```
# Check number of genes 
zgrep -c "^>" Mcavernosa.maker.transcripts.fasta
25142 

nano mcav_diamond_blastx.sh

#!/bin/bash
#SBATCH --job-name="diamond-blastx"
#SBATCH -t 240:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH -p putnamlab
#SBATCH --mem=100GB
#SBATCH --error="mcav_diamond_blastx_out_error"
#SBATCH --output="mcav_diamond_blastx_out"

echo "START" $(date)
module load DIAMOND/2.0.0-GCC-8.3.0 #Load DIAMOND

echo "Updating Mcav annotation" $(date)
diamond blastx -d /data/putnamlab/shared/databases/nr.dmnd -q Mcavernosa.maker.transcripts.fasta -o Mcav_annot -f 100 -b20 --more-sensitive -e 0.00001 -k1

echo "Search complete... converting format to XML and tab"

diamond view -a Mcav_annot.daa -o Mcav_annot.xml -f 5
diamond view -a Mcav_annot.daa -o Mcav_annot.tab -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen

echo "STOP" $(date)

sbatch mcav_diamond_blastx.sh

```

#### Ofav

```
# Check number of genes 
zgrep -c "^>" GCF_002042975.1_ofav_dov_v1_rna.fna.gz
35971 

nano ofav_diamond_blastx.sh

#!/bin/bash
#SBATCH --job-name="diamond-blastx"
#SBATCH -t 240:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH -p putnamlab
#SBATCH --mem=100GB
#SBATCH --error="ofav_diamond_blastx_out_error"
#SBATCH --output="ofav_diamond_blastx_out"

echo "START" $(date)
module load DIAMOND/2.0.0-GCC-8.3.0 #Load DIAMOND

echo "Updating Ofav annotation" $(date)
diamond blastx -d /data/putnamlab/shared/databases/nr.dmnd -q GCF_002042975.1_ofav_dov_v1_rna.fna.gz -o Ofav_annot -f 100 -b20 --more-sensitive -e 0.00001 -k1

echo "Search complete... converting format to XML and tab"

diamond view -a Ofav_annot.daa -o Ofav_annot.xml -f 5
diamond view -a Ofav_annot.daa -o Ofav_annot.tab -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen

echo "STOP" $(date)

sbatch ofav_diamond_blastx.sh

```

### Hawaii species 

#### Mcap

```
# Check number of genes 
zgrep -c "^>" Mcap.mRNA.fa
63227
```

#### Pcomp

```
# Check number of genes 
zgrep -c "^>" Porites_compressa_CDS.fa
74728
```

#### Pdam

```
# Check number of genes 
zgrep -c "^>" GCF_003704095.1_ASM370409v1_rna.fna
27287
```

#### Plob

Using plutea genome, proteins, transcripts, etc 

```
# Check number of genes 
zgrep -c "^>" plut2v1.1.transcripts.fasta
31126
```