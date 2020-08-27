### Annotating genes with R. Cunning pipeline 

MAKER pipeline [here](https://github.com/jrcunning/pdam-genome/blob/master/annotation/Makefile)

I downloaded pdam-genome repo onto uri server. Going to try to reproduce his annotation files

#### 1) Run BUSCO

Need BUSCO python script, path to pdam.fasta, and path to metazoa_odb9

```
busco/run_pdam/short_summary_pdam.txt: ../data/filter/pdam.fasta
	cd /scratch/projects/crf/pdam-genome/busco && \
	python ~/local/busco/BUSCO.py -f -c 96 --long -i ../data/filter/pdam.fasta -o pdam -l ~/local/busco/metazoa_odb9 -m geno


# obtaining busco command that Cunning used 
wget https://gitlab.com/ezlab/busco/-/raw/5b3b875377bc48f68a474813be76ed74db1f6858/BUSCO.py?inline=false
mv BUSCO.py?inline=false BUSCO.py

nano busco-run.sh

#!/bin/bash
#SBATCH --job-name="busco"
#SBATCH --time="100:00:00"
#SBATCH --nodes 1 --ntasks-per-node=20
#SBATCH --mem=250G 
##SBATCH --account=putnamlab
##SBATCH --export=NONE
#SBATCH --error="busco_out_error"
#SBATCH --output="busco_trim_out"
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu

module load BUSCO/4.0.6-foss-2019b-Python-3.7.4 

#metazoa: /data/putnamlab/shared/busco/downloads/lineages/metazoa_odb10

# Cunning didn't seem to have the fasta file on github, so i am using the one he uploaded on reef genomics. not sure it will work. as backup, use NCBI fasta file 
#fasta: /data/putnamlab/jillashey/genome/Pdam/pdam_scaffolds.fasta 


python BUSCO.py -f --long -i /data/putnamlab/jillashey/genome/Pdam/pdam_scaffolds.fasta -o /data/putnamlab/jillashey/pdam-genome/JA/busco/pdam.busco -l /data/putnamlab/shared/busco/downloads/lineages/metazoa_odb10 -m geno

sbatch busco-run.sh

```

f=force rewriting of existing files 
c=number of threads to use
long=Optimization mode Augustus self-training (Default: Off) adds considerably to the run time, but can improve results for some non-model organisms
i=input FASTA file
o=output file
l=specify what BUSCO lineage data to use
m=specify hat BUSCO analysis run to use (geno=genome assemblies)

The code above from Cunning keeps failing so I am going to use the code that Tajashree wrote for BUSCO

Got this code from /data/putnamlab/shared/busco/scripts

```
#!/bin/bash

#SBATCH --job-name="busco"
#SBATCH --time="100:00:00"
#SBATCH --nodes 1 --ntasks-per-node=20
#SBATCH --mem=250G 
##SBATCH --output="busco-%u-%x-%j"
##SBATCH --account=putnamlab
##SBATCH --export=NONE

echo "START" $(date)

labbase=/data/putnamlab
busco_shared="${labbase}/shared/busco" 
[ -z "$query" ] && query="${labbase}/jillashey/genome/Pdam/pdam_scaffolds.fasta" # set this to the query (genome/transcriptome) you are running 
[ -z "$db_to_compare" ] && db_to_compare="${busco_shared}/downloads/lineages/metazoa_odb10"

source "${busco_shared}/scripts/busco_init.sh"  # sets up the modules required for this in the right order

# we require the agustus_config/ directory copied to a "writetable" location for
# busco to run and AUGUSTUS_CONFIG_PATH set to that

if [ ! -d "${labbase}/jillashey/pdam-genome/JA/busco/agustus_config" ] ; then
    echo -e "Copying agustus_config/ to ${labbase}/jillashey/pdam-genome/JA/busco/"
    tar -C "${labbase}/jillashey/pdam-genome/JA/busco" -xzf "${busco_shared}/agustus_config.tgz"
    echo done
fi

export AUGUSTUS_CONFIG_PATH="${labbase}/${USER}/agustus_config"
# This will generate output under your $HOME/busco_output
cd "${labbase}/${USER}"
sbatch busco --config "${busco_shared}/scripts/busco-config.ini"  -f -c 20 --long -i "${query}" -l "${db_to_compare}" -o busco_output -m geno

echo "STOP" $(date)
```

#### 2) Run RepeatMasker 
