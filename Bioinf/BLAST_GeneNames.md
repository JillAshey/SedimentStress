# BLAST

Blasting genes with matching GO terms to find comparable genes from other species 

Polina supplied me with code to subset the transcript/protein sequences so I have the seqs that are associated with DEGs/GO enriched

```
module load Python/2.7.15-foss-2018b
```

## Florida

### Acerv 

FULL

1) Subset transcript/protein sequences by the significantly enriched gene names (identified through DESeq2 and GOSeq)

```
# code by Polina
import csv

proteinFile = "Acerv_assembly_v1.0.protein.fa"
proteinExport = "Acervicornis.proteins.SUBSET.fasta"
seqFile = "Acerv_assembly_v1.0.mRNA.fa"
seqExport = "Acervicornis.transcripts.SUBSET.fasta"

namesFile = "acerv_GeneNames.csv"
proteins = []
with open(namesFile) as csvfile:
    reader = csv.reader(csvfile)
    next(reader)
    for row in reader :
        proteins.append(row[1])

def scan(inFile, outFile, proteins) :
    MODE = "SEARCHING"
    with open(inFile) as PF:
        with open(outFile, 'w') as PE:
            for line in PF:
                if (line[0] == ">") :
                    heading_pname = line.split(' ')[0][1:].replace("-RA","")
                    if heading_pname in proteins :
                        MODE = "SCANNING"
                        PE.write(line)
                        proteins.remove(heading_pname)
                    else :
                        MODE = "SEARCHING"
                elif MODE == "SCANNING" :
                    PE.write(line)
    return proteins

rem_prot = scan(proteinFile, proteinExport, list(proteins))
if(len(rem_prot) > 0) :
    print("Note: Some proteins not found in protein fasta. These have been saved to protein_unfound.txt.")
    with open("protein_unfound.txt", 'w') as f:
        for item in rem_prot :
            f.write(item)

rem_seq = scan(seqFile, seqExport, list(proteins))
if(len(rem_seq) > 0) :
    print("Note: Some proteins not found in sequence fasta. These have been saved to sequence_unfound.txt.")
    with open("sequence_unfound.txt", 'w') as f:
        for item in rem_seq :
            f.write(item)
```

2) Insert a space between the > and first character of gene name 

```
sed -i -e 's/>A/> A/g' Acervicornis.transcripts.SUBSET.fasta
sed -i -e 's/>A/> A/g' Acervicornis.proteins.SUBSET.fasta
```

3) Check how many sequences in file 

```
grep -c '^>' Acervicornis.transcripts.SUBSET.fasta
# 9

grep -c '^>' Acervicornis.proteins.SUBSET.fasta
# 9 
```

4) make directory for blast output 

```
cd /data/putnamlab/jillashey/Francois_data/Florida
mkdir blast/acerv
cd blast/acerv
```

5) Run blastn

```
# nt against nt

# query fasta file
/data/putnamlab/jillashey/Francois_data/Florida/blast/acerv/Acervicornis.transcripts.SUBSET.fasta

# Will blast remotely - NCBI servers

nano acerv_transcript_blast.sh

#!/bin/bash
#SBATCH --job-name="blastn"
#SBATCH -t 30-00:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH -D /data/putnamlab/jillashey/Francois_data/Florida/blast/acerv
#SBATCH --error="acerv_transcript_remote_blast_out_error"
#SBATCH --output="acerv_transcript_remote_blast_out"

module load BLAST+/2.8.1-foss-2018b 

blastn -db nt -query Acervicornis.transcripts.SUBSET.fasta -out acerv_transcript_blast.out -remote

sbatch acerv_transcript_blast.sh
Submitted batch job 1837022
```

Running blastn on all acerv transcript seqs

```
# nt against nt

# query fasta file
/data/putnamlab/jillashey/Francois_data/Florida/blast/acerv/Acervicornis.transcripts.SUBSET.fasta

# Will blast remotely - NCBI servers

# Running on putnam lab node, 

nano acerv_AllTranscript_blast.sh

#!/bin/bash
#SBATCH --job-name="blastn"
#SBATCH -t 30-00:00:00
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH -D /data/putnamlab/jillashey/Francois_data/Florida/blast/acerv
#SBATCH --error="acerv_AllTranscript_remote_blast_out_error"
#SBATCH --output="acerv_AllTranscript_remote_blast_out"

module load BLAST+/2.8.1-foss-2018b

blastn -db nt -query Acerv_assembly_v1.0.mRNA.fa -out acerv_AllTranscript_blast.out -max_target_seqs 2 -outfmt 7 -remote 

sbatch acerv_AllTranscript_blast.sh
Submitted batch job 1846442
```

6) Run blastp

```
# prot against prot

# query fasta file
/data/putnamlab/jillashey/Francois_data/Florida/blast/acerv/Acervicornis.proteins.SUBSET.fasta

# Will blast remotely - NCBI servers

nano acerv_protein_blast.sh

#!/bin/bash
#SBATCH --job-name="blastp"
#SBATCH -t 30-00:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH -D /data/putnamlab/jillashey/Francois_data/Florida/blast/acerv
#SBATCH --error="acerv_protein_remote_blast_out_error"
#SBATCH --output="acerv_protein_remote_blast_out"

module load BLAST+/2.8.1-foss-2018b 

blastp -db nr -query Acervicornis.proteins.SUBSET.fasta -remote -out acerv_protein_blast.out 

sbatch acerv_protein_blast.sh
Submitted batch job 1837023

```

SUBSET of Acerv - removed outlier samples 

1) Subset transcript/protein sequences by the significantly enriched gene names (identified through DESeq2 and GOSeq)

```
# code by Polina
import csv

proteinFile = "Acerv_assembly_v1.0.protein.fa"
proteinExport = "Acervicornis_sub.proteins.SUBSET.fasta"
seqFile = "Acerv_assembly_v1.0.mRNA.fa"
seqExport = "Acervicornis_sub.transcripts.SUBSET.fasta"

namesFile = "acerv_sub_GeneNames.csv"
proteins = []
with open(namesFile) as csvfile:
    reader = csv.reader(csvfile)
    next(reader)
    for row in reader :
        proteins.append(row[1])

def scan(inFile, outFile, proteins) :
    MODE = "SEARCHING"
    with open(inFile) as PF:
        with open(outFile, 'w') as PE:
            for line in PF:
                if (line[0] == ">") :
                    heading_pname = line.split(' ')[0][1:].replace("-RA","")
                    if heading_pname in proteins :
                        MODE = "SCANNING"
                        PE.write(line)
                        proteins.remove(heading_pname)
                    else :
                        MODE = "SEARCHING"
                elif MODE == "SCANNING" :
                    PE.write(line)
    return proteins

rem_prot = scan(proteinFile, proteinExport, list(proteins))
if(len(rem_prot) > 0) :
    print("Note: Some proteins not found in protein fasta. These have been saved to protein_unfound.txt.")
    with open("protein_unfound.txt", 'w') as f:
        for item in rem_prot :
            f.write(item)

rem_seq = scan(seqFile, seqExport, list(proteins))
if(len(rem_seq) > 0) :
    print("Note: Some proteins not found in sequence fasta. These have been saved to sequence_unfound.txt.")
    with open("sequence_unfound.txt", 'w') as f:
        for item in rem_seq :
            f.write(item)
```

2) Insert a space between the > and first character of gene name 

```
sed -i -e 's/>A/> A/g' Acervicornis_sub.transcripts.SUBSET.fasta
sed -i -e 's/>A/> A/g' Acervicornis_sub.proteins.SUBSET.fasta
```

3) Check how many sequences in file 

```
grep -c '^>' Acervicornis_sub.transcripts.SUBSET.fasta
# 39

grep -c '^>' Acervicornis_sub.transcripts.SUBSET.fasta
# 39
```

4) make directory for blast output 

Acerv blast folder made above

5) Run blastn

```
# nt against nt

# query fasta file
/data/putnamlab/jillashey/Francois_data/Florida/blast/acerv/Acervicornis.transcripts.SUBSET.fasta

# Will blast remotely - NCBI servers

nano acerv_transcript_sub_blast.sh

#!/bin/bash
#SBATCH --job-name="acerv_blastn"
#SBATCH -t 100-00:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH -D /data/putnamlab/jillashey/Francois_data/Florida/blast/acerv
#SBATCH --error="acerv_transcript_sub_remote_blast_out_error"
#SBATCH --output="acerv_transcript_sub_remote_blast_out"

module load BLAST+/2.8.1-foss-2018b 

blastn -db nt -query Acervicornis_sub.transcripts.SUBSET.fasta -out acerv_transcript_sub_blast.out -max_target_seqs 2 -remote

sbatch acerv_transcript_sub_blast.sh
Submitted batch job 1845864
```



### Mcav

1) Subset transcript/protein sequences by the significantly enriched gene names (identified through DESeq2 and GOSeq)

```
import csv

proteinFile = "Mcavernosa.maker.proteins.fasta"
proteinExport = "Mcavernosa.maker.proteins.SUBSET.fasta"
seqFile = "Mcavernosa.maker.transcripts.fasta"
seqExport = "Mcavernosa.maker.transcripts.SUBSET.fasta"

namesFile = "mcav_GeneNames.csv"
proteins = []
with open(namesFile) as csvfile:
    reader = csv.reader(csvfile)
    next(reader)
    for row in reader :
        proteins.append(row[1])

def scan(inFile, outFile, proteins) :
    MODE = "SEARCHING"
    with open(inFile) as PF:
        with open(outFile, 'w') as PE:
            for line in PF:
                if (line[0] == ">") :
                    heading_pname = line.split(' ')[0][1:].replace("-RA","")
                    if heading_pname in proteins :
                        MODE = "SCANNING"
                        PE.write(line)
                        proteins.remove(heading_pname)
                    else :
                        MODE = "SEARCHING"
                elif MODE == "SCANNING" :
                    PE.write(line)
    return proteins

rem_prot = scan(proteinFile, proteinExport, list(proteins))
if(len(rem_prot) > 0) :
    print("Note: Some proteins not found in protein fasta. These have been saved to protein_unfound.txt.")
    with open("protein_unfound.txt", 'w') as f:
        for item in rem_prot :
            f.write(item)

rem_seq = scan(seqFile, seqExport, list(proteins))
if(len(rem_seq) > 0) :
    print("Note: Some proteins not found in sequence fasta. These have been saved to sequence_unfound.txt.")
    with open("sequence_unfound.txt", 'w') as f:
        for item in rem_seq :
            f.write(item)
```


2) Insert a space between the > and first character of gene name 

```
sed -i -e 's/>M/> M/g' Mcavernosa.maker.transcripts.SUBSET.fasta
sed -i -e 's/>M/> M/g' Mcavernosa.maker.proteins.SUBSET.fasta

```

3) Check how many sequences in file 

```
grep -c '^>' Mcavernosa.maker.transcripts.SUBSET.fasta
# 22

grep -c '^>' Mcavernosa.maker.proteins.SUBSET.fasta
# 22
```


4) Make directory for blast output 

```
cd /data/putnamlab/jillashey/Francois_data/Florida/blast
mkdir mcav
cd mcav
```

5) Run blastn

```
# nt against nt

# query fasta file
/data/putnamlab/jillashey/Francois_data/Florida/blast/mcav/Mcavernosa.maker.transcripts.SUBSET.fasta

# Will blast remotely - NCBI servers

nano mcav_transcript_blast.sh

#!/bin/bash
#SBATCH --job-name="MC_blastn"
#SBATCH -t 30-00:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH -D /data/putnamlab/jillashey/Francois_data/Florida/blast/mcav
#SBATCH --error="mcav_transcript_remote_blast_out_error"
#SBATCH --output="mcav_transcript_remote_blast_out"

module load BLAST+/2.8.1-foss-2018b 

blastn -db nt -query Mcavernosa.maker.transcripts.SUBSET.fasta -out mcav_transcript_blast.out -remote

sbatch mcav_transcript_blast.sh
Submitted batch job 1845863

```

Editing script so that blast only returns 2 target seqs max for each query 

```
# nt against nt

# query fasta file
/data/putnamlab/jillashey/Francois_data/Florida/blast/mcav/Mcavernosa.maker.transcripts.SUBSET.fasta

# Will blast remotely - NCBI servers

nano mcav_transcript_blast_maxseqs.sh

#!/bin/bash
#SBATCH --job-name="blastn"
#SBATCH -t 30-00:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH -D /data/putnamlab/jillashey/Francois_data/Florida/blast/mcav
#SBATCH --error="mcav_transcript_remote_blast_maxseqs_out_error"
#SBATCH --output="mcav_transcript_remote_blast_maxseqs_out"

module load BLAST+/2.8.1-foss-2018b 

blastn -db nt -query Mcavernosa.maker.transcripts.SUBSET.fasta -out mcav_transcript_maxseqs_blast.out -max_target_seqs 2 -remote 

sbatch mcav_transcript_blast_maxseqs.sh
Submitted batch job 1845863

```

Running blastn on all mcav transcript seqs

```
# nt against nt

# query fasta file
/data/putnamlab/jillashey/Francois_data/Florida/blast/mcav/Mcavernosa.maker.transcripts.fasta

# Will blast remotely - NCBI servers

# Running on putnam lab node, 

nano mcav_AllTranscript_blast.sh

#!/bin/bash
#SBATCH --job-name="blastn"
#SBATCH -t 100-00:00:00
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH -D /data/putnamlab/jillashey/Francois_data/Florida/blast/mcav
#SBATCH --error="mcav_AllTranscript_remote_blast_out_error"
#SBATCH --output="mcav_AllTranscript_remote_blast_out"

module load BLAST+/2.8.1-foss-2018b

blastn -db nt -query Mcavernosa.maker.transcripts.fasta -out mcav_AllTranscript_blast.out -max_target_seqs 2 -outfmt 7 -remote 

sbatch mcav_AllTranscript_blast.sh
Submitted batch job 1846444
```


6) Run blastp

```
# prot against prot

# query fasta file
/data/putnamlab/jillashey/Francois_data/Florida/blast/mcav/Mcavernosa.maker.proteins.SUBSET.fasta

# Will blast remotely - NCBI servers

nano mcav_protein_blast.sh

#!/bin/bash
#SBATCH --job-name="blastp"
#SBATCH -t 30-00:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH -D /data/putnamlab/jillashey/Francois_data/Florida/blast/
#SBATCH --error="mcav_protein_remote_blast_out_error"
#SBATCH --output="mcav_protein_remote_blast_out"

module load BLAST+/2.8.1-foss-2018b 

blastp -db nr -query Mcavernosa.maker.proteins.SUBSET.fasta -remote -out mcav_protein_blast.out 

sbatch mcav_protein_blast.sh
Submitted batch job 1836971

```

## Hawaii 

### Pdam 

```
import csv

proteinFile = "GCF_003704095.1_ASM370409v1_protein.faa"
proteinExport = "Pdam.proteins.SUBSET.fasta"
seqFile = "GCF_003704095.1_ASM370409v1_rna.fna"
seqExport = "Pdam.transcripts.SUBSET.fasta"

namesFile = "pdam_GeneNames.csv"
proteins = []
with open(namesFile) as csvfile:
    reader = csv.reader(csvfile)
    next(reader)
    for row in reader :
        proteins.append(row[1])

def scan(inFile, outFile, proteins) :
    MODE = "SEARCHING"
    with open(inFile) as PF:
        with open(outFile, 'w') as PE:
            for line in PF:
                if (line[0] == ">") :
                    heading_pname = line.split(' ')[0][1:].replace("-RA","")
                    if heading_pname in proteins :
                        MODE = "SCANNING"
                        PE.write(line)
                        proteins.remove(heading_pname)
                    else :
                        MODE = "SEARCHING"
                elif MODE == "SCANNING" :
                    PE.write(line)
    return proteins

rem_prot = scan(proteinFile, proteinExport, list(proteins))
if(len(rem_prot) > 0) :
    print("Note: Some proteins not found in protein fasta. These have been saved to protein_unfound.txt.")
    with open("protein_unfound.txt", 'w') as f:
        for item in rem_prot :
            f.write(item)

rem_seq = scan(seqFile, seqExport, list(proteins))
if(len(rem_seq) > 0) :
    print("Note: Some proteins not found in sequence fasta. These have been saved to sequence_unfound.txt.")
    with open("sequence_unfound.txt", 'w') as f:
        for item in rem_seq :
            f.write(item)
```




Running blastn on all pdam transcript seqs

##### NCBI

```
# nt against nt

# query fasta file
ln -s /data/putnamlab/jillashey/genome/Pdam/NCBI/GCF_003704095.1_ASM370409v1_rna.fna

# Will blast remotely - NCBI servers

nano pdam_AllTranscript_blast.sh

#!/bin/bash
#SBATCH --job-name="blastn"
#SBATCH -t 100-00:00:00
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH -D /data/putnamlab/jillashey/Francois_data/Hawaii/blast/pdam/NCBI
#SBATCH --error="pdam_AllTranscript_remote_blast_out_error"
#SBATCH --output="pdam_AllTranscript_remote_blast_out"

module load BLAST+/2.8.1-foss-2018b

blastn -db nt -query GCF_003704095.1_ASM370409v1_rna.fna -out pdam_AllTranscript_blast.out -max_target_seqs 2 -outfmt 7 -remote 

sbatch pdam_AllTranscript_blast.sh
Submitted batch job 1846514
```

##### Reef Genomics

```
# nt against nt

# query fasta file
ln -s /data/putnamlab/jillashey/genome/Pdam/ReefGenomics/pdam_transcripts.fasta

# Will blast remotely - NCBI servers

nano pdam_RG_AllTranscript_blast.sh

#!/bin/bash
#SBATCH --job-name="blastn"
#SBATCH -t 100-00:00:00
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH -D /data/putnamlab/jillashey/Francois_data/Hawaii/blast/pdam/ReefGenomics
#SBATCH --error="pdam_RG_AllTranscript_remote_blast_out_error"
#SBATCH --output="pdam_RG_AllTranscript_remote_blast_out"

module load BLAST+/2.8.1-foss-2018b

blastn -db nt -query pdam_transcripts.fasta -out pdam_RG_AllTranscript_blast.out -max_target_seqs 2 -outfmt 7 -remote 

sbatch pdam_RG_AllTranscript_blast.sh
Submitted batch job 1846515
```



### Plob

1) Subset transcript/protein sequences by the significantly enriched gene names (identified through DESeq2 and GOSeq)

```
import csv

proteinFile = "plut2v1.1.proteins.fasta"
proteinExport = "Plobata.proteins.SUBSET.fasta"
seqFile = "plut2v1.1.transcripts.fasta"
seqExport = "Plobata.transcripts.SUBSET.fasta"

namesFile = "plob_GeneNames.csv"
proteins = []
with open(namesFile) as csvfile:
    reader = csv.reader(csvfile)
    next(reader)
    for row in reader :
        proteins.append(row[1])

def scan(inFile, outFile, proteins) :
    MODE = "SEARCHING"
    with open(inFile) as PF:
        with open(outFile, 'w') as PE:
            for line in PF:
                if (line[0] == ">") :
                    heading_pname = line.split(' ')[0][1:].replace("-RA","")
                    if heading_pname in proteins :
                        MODE = "SCANNING"
                        PE.write(line)
                        proteins.remove(heading_pname)
                    else :
                        MODE = "SEARCHING"
                elif MODE == "SCANNING" :
                    PE.write(line)
    return proteins

rem_prot = scan(proteinFile, proteinExport, list(proteins))
if(len(rem_prot) > 0) :
    print("Note: Some proteins not found in protein fasta. These have been saved to protein_unfound.txt.")
    with open("protein_unfound.txt", 'w') as f:
        for item in rem_prot :
            f.write(item)

rem_seq = scan(seqFile, seqExport, list(proteins))
if(len(rem_seq) > 0) :
    print("Note: Some proteins not found in sequence fasta. These have been saved to sequence_unfound.txt.")
    with open("sequence_unfound.txt", 'w') as f:
        for item in rem_seq :
            f.write(item)
```




Running blastn on all plob transcript seqs

```
# nt against nt

# query fasta file
plut2v1.1.transcripts.fasta


ln -s /data/putnamlab/jillashey/genome/Plutea/plut2v1.1.transcripts.fasta

# Will blast remotely - NCBI servers

nano plob_AllTranscript_blast.sh

#!/bin/bash
#SBATCH --job-name="blastn"
#SBATCH -t 100-00:00:00
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH -D /data/putnamlab/jillashey/Francois_data/Hawaii/blast/plob
#SBATCH --error="plob_AllTranscript_remote_blast_out_error"
#SBATCH --output="plob_AllTranscript_remote_blast_out"

module load BLAST+/2.8.1-foss-2018b

blastn -db nt -query plut2v1.1.transcripts.fasta -out plob_AllTranscript_blast.out -max_target_seqs 2 -outfmt 7 -remote 

sbatch plob_AllTranscript_blast.sh
Submitted batch job 1846518
```

