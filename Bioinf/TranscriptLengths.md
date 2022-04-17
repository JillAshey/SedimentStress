In bluewaves, calculate gene lengths 

## Florida

### Acerv 

```
# Transcript lengths
awk 'BEGIN{FS="[> ]"} /^>/{val=$2;next}  {print val,length($0)}' Acerv_assembly_v1.0.mRNA.fa > length.mRNA_Acerv.txt

# Protein lengths 
awk 'BEGIN{FS="[> ]"} /^>/{val=$2;next}  {print val,length($0)}' Acerv_assembly_v1.0.protein.fa > length.protein_Acerv.txt

```

### Mcav 

Sequences were tab-deliminated in the fasta file. I had to convert the 'multi-line' fasta to a single line fasta 

```
awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' Mcavernosa.maker.transcripts.fasta > newfile.fasta

awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' Mcavernosa.maker.proteins.SUBSET.fasta > mcav_newfile_protein.fasta

```

Now gene lengths can be calculated 

```
awk 'BEGIN{FS="[> ]"} /^>/{val=$2;next}  {print val,length($0)}' newfile.fasta > length.mRNA_Acerv.txt

```

### Ofav 

Sequences were tab-deliminated in the fasta file. I had to convert the 'multi-line' fasta to a single line fasta 

```
awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' GCF_002042975.1_ofav_dov_v1_rna.fna > newfile_RNA_Ofav.fasta

awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' GCF_002042975.1_ofav_dov_v1_protein.faa > newfile_Protein_Ofav.fasta
```

Now gene lengths can be calculated 

```
# Transcript lengths
awk 'BEGIN{FS="[>]"} /^>/{val=$2;next}  {print val,length($0)}' newfile_RNA_Ofav.fasta > length.mRNA_Ofav.txt

# Protein lengths 
awk 'BEGIN{FS="[>]"} /^>/{val=$2;next}  {print val,length($0)}' newfile_Protein_Ofav.fasta > length.protein_Ofav.txt

```

## Hawaii

### Mcap

Sequences were tab-deliminated in the fasta file. I had to convert the 'multi-line' fasta to a single line fasta 

```
awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' Montipora_capitata_HIv1.genesNoCopies.gff3_polished.cds.fna > newfile_RNA_Mcap.fasta

awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' Montipora_capitata_HIv1.genesNoCopies.gff3_polished.pep.faa > newfile_Protein_Mcap.fasta
```

Now gene lengths can be calculated 

```
# Transcript lengths
awk 'BEGIN{FS="[> ]"} /^>/{val=$2;next}  {print val,length($0)}' newfile_RNA_Mcap.fasta > length.mRNA_Mcap.txt

# Protein lengths 
awk 'BEGIN{FS="[> ]"} /^>/{val=$2;next}  {print val,length($0)}' newfile_Protein_Mcap.fasta > length.protein_Mcap.txt
```

### Pacuta 

Sequences were tab-deliminated in the fasta file. I had to convert the 'multi-line' fasta to a single line fasta 

```
awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' braker_v1.codingseq.fasta > newfile_RNA_Pacuta.fasta

awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' Pocillopora_acuta_HIv1.genes.pep.faa > newfile_Protein_Pacuta.fasta

```

Now gene lengths can be calculated 

```
# Transcript lengths
awk 'BEGIN{FS="[> ]"} /^>/{val=$2;next}  {print val,length($0)}' newfile_RNA_Pacuta.fasta > length.mRNA_Pacuta.txt

# Protein lengths 
awk 'BEGIN{FS="[> ]"} /^>/{val=$2;next}  {print val,length($0)}' newfile_Protein_Pacuta.fasta > length.protein_Pacuta.txt
```

### Plob 

Sequences were tab-deliminated in the fasta file. I had to convert the 'multi-line' fasta to a single line fasta 

```
awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' plut2v1.1.transcripts.fasta > newfile_mRNA_Plob.fasta

awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' plut2v1.1.proteins.fasta > newfile_Protein_Plob.fasta

```

Now gene lengths can be calculated 

```
# Transcript lengths
awk 'BEGIN{FS="[> ]"} /^>/{val=$2;next}  {print val,length($0)}' newfile_mRNA_Plob.fasta > length.mRNA_Plob.txt

# Protein lengths 
awk 'BEGIN{FS="[> ]"} /^>/{val=$2;next}  {print val,length($0)}' newfile_Protein_Plob.fasta > length.protein_Plob.txt
```