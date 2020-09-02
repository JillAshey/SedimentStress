### Annotating genes with R. Cunning pipeline 

MAKER pipeline [here](https://github.com/jrcunning/pdam-genome/blob/master/annotation/Makefile)

I downloaded pdam-genome repo onto uri server. Going to try to reproduce MAKER annotation files. I will use the files from Reef Genomics to begin.

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

Yay BUSCO ran successfully! 

```
# BUSCO version is: 4.0.6 
# The lineage dataset is: metazoa_odb10 (Creation date: 2019-11-20, number of species: 65, number of BUSCOs: 954)
# Summarized benchmarking in BUSCO notation for file /data/putnamlab/jillashey/genome/Pdam/pdam_scaffolds.fasta
# BUSCO was run in mode: genome

        ***** Results: *****

        C:89.3%[S:88.9%,D:0.4%],F:3.4%,M:7.3%,n:954        
        852     Complete BUSCOs (C)                        
        848     Complete and single-copy BUSCOs (S)        
        4       Complete and duplicated BUSCOs (D)         
        32      Fragmented BUSCOs (F)                      
        70      Missing BUSCOs (M)                         
        954     Total BUSCO groups searched  
```

#### 2) Obtain data for maker

The next step in Cunning pipeline is to get all data in order for maker. Looking at his maker ctl files, I don't have the same files at all to put into the maker control files. I can probably generate these myself.

Files I need to make/obtain: 

- evidence from two RNAseq datasets from Pdam (he used (A)Traylor-Knowles, N. et al. Production of a reference transcriptome and transcriptomic database (PocilloporaBase) for the cauliflower coral, Pocillopora damicornis. BMC Genomics 12, 585 (2011). (B) Mayfield, A. B., Wang, Y.-B., Chen, C.-S., Lin, C.-Y. & Chen, S.-H. Compartment-specific transcriptomics in a reef-building coral exposed)
	- I didn't obtain his exact files, but pretty close I think. On server under /data/putnamlab/jillashey/genome/Pdam/ReefGenomics/pdam_transcripts.fasta and Pocillopora_damicornis_cds_100.final.clstr.fna
- Evidence from closely related species - S. pistillata (Voolstra, C. R. et al. Comparative analysis of the genomes of Stylophora pistillata and Acropora digitifera provides evidence for extensive differences between species of corals. Sci. Rep. 7, 17583 (2017).)
	- Found it! Using data from Reef Genomics. On server, under /data/putnamlab/jillashey/genome/Spist/Spis.Trinity.cov10.longest.fa.gz
- Protein seqs from 20 coral species (Bhattacharya, D. et al. Comparative genomics explains the evolutionary success of reef-forming corals. Elife 5 (2016).)
	- Lots of species info on Reef Genomics. Don't see a file containing all coral protein seqs. Not sure it would work to combine them all into one file, but I could also just list them out
	- Using species on server (in /data/putnamlab/jillashey/genome/)
		- Adig (Adig/GCF_000222465.1_Adig_1.1_protein.faa.gz)
		- Amill (Amil_v2.01/Amil.all.maker.proteins.fasta ?)
		- Aten (Aten/aten_0.11.maker_post_001.proteins.fasta.gz)
		- Gfasc (Gfasc/gfas_1.0.proteins.fasta.gz)
		- Mcap (Mcap/Mcap.protein.fa.gz)
		- Pdam (Pdam/ReefGenomics/pdam_proteins.fasta.gz)
		- Spist (Spist/Spis.genome.annotation.pep.longest.fa.gz)
- rmLib - organism specific repeat library in fasta format 
	- Ross had this in his github. On server under /data/putnamlab/jillashey/pdam-genome/annotation/repeatmasker/consensi.fa.classified
- repeat_protein - te_proteins.fasta - file of transposable element proteins 
- augustus_species - BUSCO_pdam_1629150601 - augustus gene prediction model, likely with BUSCO output

To get te_proteins, I think I need to run RepeatModeler

#### 3) Run RepeatModeler

My ```RepeatModeler``` script wasn't working, so I downloaded RepeatModeler myself.

/data/putnamlab/jillashey/pdam-genome/JA/RepeatModeler/RepeatModeler-2.0.1/

```
wget http://www.repeatmasker.org/RepeatModeler/RepeatModeler-2.0.1.tar.gz

tar -zxvf RepeatModeler-2.0.1.tar.gz

cd RepeatModeler-2.0.1

```

After you actually get into RepeatModeler, it will ask you to configure a bunch of different modules

```
perl configure 

# Will prompt you to specify path to various software packages
perl:/opt/software/Perl/5.24.0-foss-2016b/bin/perl
RepeatMasker:/opt/software/RepeatMasker/
RECON de-novo repeatfinding program:/opt/RECON-1.08/bin
RepeatScout de-novo repeatfinding program (1.0.6 or higher): /opt/RepeatScout-1.0.5 # was under the specified version #
TRF program ( 4.0.9 or higher): /opt/software/TRF/4.09-linux64/

easybuild-TRF-4.09-20170106.225221.log  easybuild-TRF-4.09-20170106.225221_test_report.md  TRF-4.09-linux64-easybuild-devel  TRF-4.09-linux64.eb

trf

```

May just skip repeat masking step for now

#### 4) Run maker 

```
module load maker/3.01.03
maker -CTL # generates files to use with maker 
```

I'm going to do a test run of maker before doing Ross code. Also need to find some parts to truly replicate his code 

```
#-----Genome (these are always required)
genome=data/putnamlab/jillashey/genome/Pdam/pdam_scaffolds.fasta #genome sequence (fasta file or fasta embeded in GFF3 file)
organism_type=eukaryotic #eukaryotic or prokaryotic. Default is eukaryotic

#-----Re-annotation Using MAKER Derived GFF3
maker_gff= #MAKER derived GFF3 file
est_pass=0 #use ESTs in maker_gff: 1 = yes, 0 = no
altest_pass=0 #use alternate organism ESTs in maker_gff: 1 = yes, 0 = no
protein_pass=0 #use protein alignments in maker_gff: 1 = yes, 0 = no
rm_pass=0 #use repeats in maker_gff: 1 = yes, 0 = no
model_pass=0 #use gene models in maker_gff: 1 = yes, 0 = no
pred_pass=0 #use ab-initio predictions in maker_gff: 1 = yes, 0 = no
other_pass=0 #passthrough anyything else in maker_gff: 1 = yes, 0 = no

#-----EST Evidence (for best results provide a file for at least one)
est=data/putnamlab/jillashey/genome/Pdam/ReefGenomics/pdam_transcripts.fasta,data/putnamlab/jillashey/genome/Pdam/ReefGenomics/Pocillopora_damicornis_cds_100.final.clstr.fna #set of ESTs or assembled mRNA-seq in fasta format
altest=data/putnamlab/jillashey/genome/Spist/Spis.Trinity.cov10.longest.fa #EST/cDNA sequence file in fasta format from an alternate organism
est_gff= #aligned ESTs or mRNA-seq from an external GFF3 file
altest_gff= #aligned ESTs from a closly relate species in GFF3 format

#-----Protein Homology Evidence (for best results provide a file for at least one)
protein=data/putnamlab/jillashey/genome/Pdam/ReefGenomics/pdam_proteins.fasta,data/putnamlab/jillashey/genome/Acerv/Acerv_assembly_v1.0.protein.fa,data/putnamlab/jillashey/genome/Adig/GCF_000222465.1_Adig_1.1_protein.faa,data/putnamlab/jillashey/genome/Amil_v2.01/Amil.all.maker.proteins.fasta,data/putnamlab/jillashey/genome/Aten/aten_0.11.maker_post_001.proteins.fasta,data/putnamlab/jillashey/genome/Gfasc/gfas_1.0.proteins.fasta,data/putnamlab/jillashey/genome/Mcap/Mcap.protein.fa,data/putnamlab/jillashey/genome/Mcav/Mcav_genome/Mcavernosa_annotation/Mcavernosa.maker.proteins.fasta #protein sequence file in fasta format (i.e. from mutiple organisms)
protein_gff=  #aligned protein homology evidence from an external GFF3 file

#-----Repeat Masking (leave values blank to skip repeat masking)
model_org=all #select a model organism for RepBase masking in RepeatMasker
rmlib=data/putnamlab/jillashey/pdam-genome/annotation/repeatmasker/consensi.fa.classified #provide an organism specific repeat library in fasta format for RepeatMasker
repeat_protein= #provide a fasta file of transposable element proteins for RepeatRunner
rm_gff= #pre-identified repeat elements from an external GFF3 file
prok_rm=0 #forces MAKER to repeatmask prokaryotes (no reason to change this), 1 = yes, 0 = no
softmask=1 #use soft-masking rather than hard-masking in BLAST (i.e. seg and dust filtering)

#-----Gene Prediction
snaphmm= #SNAP HMM file
gmhmm= #GeneMark HMM file
augustus_species= #Augustus gene prediction species model
fgenesh_par_file= #FGENESH parameter file
pred_gff= #ab-initio predictions from an external GFF3 file
model_gff= #annotated gene models from an external GFF3 file (annotation pass-through)
run_evm=0 #run EvidenceModeler, 1 = yes, 0 = no
est2genome=0 #infer gene predictions directly from ESTs, 1 = yes, 0 = no
protein2genome=0 #infer predictions from protein homology, 1 = yes, 0 = no
trna=0 #find tRNAs with tRNAscan, 1 = yes, 0 = no
snoscan_rrna= #rRNA file to have Snoscan find snoRNAs
snoscan_meth= #-O-methylation site fileto have Snoscan find snoRNAs
unmask=0 #also run ab-initio prediction programs on unmasked sequence, 1 = yes, 0 = no
allow_overlap= #allowed gene overlap fraction (value from 0 to 1, blank for default)

#-----Other Annotation Feature Types (features MAKER doesn't recognize)
other_gff= #extra features to pass-through to final MAKER generated GFF3 file

#-----External Application Behavior Options
alt_peptide=C #amino acid used to replace non-standard amino acids in BLAST databases
cpus=1 #max number of cpus to use in BLAST and RepeatMasker (not for MPI, leave 1 when using MPI)

#-----MAKER Behavior Options
max_dna_len=100000 #length for dividing up contigs into chunks (increases/decreases memory usage)
min_contig=1 #skip genome contigs below this length (under 10kb are often useless)

pred_flank=200 #flank for extending evidence clusters sent to gene predictors
pred_stats=0 #report AED and QI statistics for all predictions as well as models
AED_threshold=1 #Maximum Annotation Edit Distance allowed (bound by 0 and 1)
min_protein=0 #require at least this many amino acids in predicted proteins
alt_splice=0 #Take extra steps to try and find alternative splicing, 1 = yes, 0 = no
always_complete=0 #extra steps to force start and stop codons, 1 = yes, 0 = no
map_forward=0 #map names and attributes forward from old GFF3 genes, 1 = yes, 0 = no
keep_preds=0 #Concordance threshold to add unsupported gene prediction (bound by 0 and 1)

split_hit=10000 #length for the splitting of hits (expected max intron size for evidence alignments)
min_intron=20 #minimum intron length (used for alignment polishing)
single_exon=0 #consider single exon EST evidence when generating annotations, 1 = yes, 0 = no
single_length=250 #min length required for single exon ESTs if 'single_exon is enabled'
correct_est_fusion=0 #limits use of ESTs in annotation to avoid fusion genes

tries=2 #number of times to try a contig if there is a failure for some reason
clean_try=0 #remove all data from previous run before retrying, 1 = yes, 0 = no
clean_up=0 #removes theVoid directory with individual analysis files, 1 = yes, 0 = no
TMP= #specify a directory other than the system default temporary directory for temporary files
```

```
nano maker_round0.sh

#!/bin/bash
#SBATCH --job-name="maker"
#SBATCH --time="100:00:00"
#SBATCH --nodes 1 --ntasks-per-node=20
#SBATCH --mem=250G 
##SBATCH --account=putnamlab
##SBATCH -D /data/putnamlab/jillashey/pdam-genome/JA/maker
##SBATCH --export=NONE
#SBATCH --error="maker0_out_error"
#SBATCH --output="maker0_trim_out"
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu

module load maker/3.01.03


mpiexec -n 20 maker -base pdam maker_opts.ctl maker_bopts.ctl maker_exe.ctl -fix_nucleotides

sbatch maker_round0.sh

# ended up running it like this instead of .sh
sbatch mpiexec -n 20 maker -base pdam maker_opts.ctl maker_bopts.ctl maker_exe.ctl -fix_nucleotides

```

Got this error: 
Traceback (most recent call last):
  File "/var/spool/slurmd/job1689772/slurm_script", line 186, in <module>
    import mpdlib
ImportError: No module named mpdlib

```
# Attempting to install MAKER on my own 
wget http://yandell.topaz.genetics.utah.edu/maker_downloads/5E72/7851/9BF9/4C946835C1B78781970D33DE04E5/maker-3.01.03.tgz

gunzip maker-3.01.03.tgz
tar -xvf maker-3.01.03.tar
```


Use INSTALL file for installation instructions
```
cd src
perl Build.PL

Checking prerequisites...
  recommends:
    *  DBD::Pg is not installed

ERRORS/WARNINGS FOUND IN PREREQUISITES.  You may wish to install the versions
of the modules indicated above before proceeding with this installation

Run 'Build installdeps' to install missing prerequisites.


MAKER supports distributed parallelization via MPI.
Would you like to configure MAKER for MPI (This
requires that you have an MPI client installed)? [N ]
N 
Created MYMETA.yml and MYMETA.json
Creating new 'Build' script for 'MAKER' version 'v3.1.3'


The file 'Build' has been created for you to finish installing MAKER.


==============================================================================
STATUS MAKER v3.1.3
==============================================================================
PERL Dependencies:	VERIFIED
External Programs:	VERIFIED
External C Libraries:	VERIFIED
MPI SUPPORT:		DISABLED
MWAS Web Interface:	DISABLED
MAKER PACKAGE:		CONFIGURATION OK


Important Commands:
	./Build installdeps	#installs missing PERL dependencies
	./Build installexes	#installs all missing external programs
	./Build install		#installs MAKER
	./Build status		#Shows this status menu

Other Commands:
	./Build repeatmasker	#installs RepeatMasker (asks for RepBase)
	./Build blast		#installs BLAST (NCBI BLAST+)
	./Build exonerate	#installs Exonerate (v2 on UNIX / v1 on Mac OSX)
	./Build snap		#installs SNAP
	./Build augustus	#installs Augustus
	./Build apollo		#installs Apollo
	./Build gbrowse		#installs GBrowse (must be root)
	./Build jbrowse		#installs JBrowse (MAKER copy, not web accecible)
	./Build webapollo	#installs WebApollo (use maker2wap to create DBs)
	./Build mpich2		#installs MPICH2 (but manual install recommended)

./Build install
Building MAKER
Installing MAKER...
Building MAKER
Installing /data/putnamlab/jillashey/pdam-genome/JA/maker/src/../perl/lib/MAKER/ConfigData.pm
Installing /data/putnamlab/jillashey/pdam-genome/JA/maker/src/../perl/lib/Parallel/Application/MPI.pm
Installing /data/putnamlab/jillashey/pdam-genome/JA/maker/src/../perl/man/MAKER::ConfigData.3
Installing /data/putnamlab/jillashey/pdam-genome/JA/maker/src/../bin/compare_gff3_to_chado
Installing /data/putnamlab/jillashey/pdam-genome/JA/maker/src/../bin/cufflinks2gff3
Installing /data/putnamlab/jillashey/pdam-genome/JA/maker/src/../bin/evaluator
Installing /data/putnamlab/jillashey/pdam-genome/JA/maker/src/../bin/iprscan2gff3
Installing /data/putnamlab/jillashey/pdam-genome/JA/maker/src/../bin/maker2jbrowse
Installing /data/putnamlab/jillashey/pdam-genome/JA/maker/src/../bin/maker_functional_fasta
Installing /data/putnamlab/jillashey/pdam-genome/JA/maker/src/../bin/maker_functional_gff
Installing /data/putnamlab/jillashey/pdam-genome/JA/maker/src/../bin/map_fasta_ids
Installing /data/putnamlab/jillashey/pdam-genome/JA/maker/src/../bin/match2gene.pl
Installing /data/putnamlab/jillashey/pdam-genome/JA/maker/src/../bin/quality_filter.pl
Installing /data/putnamlab/jillashey/pdam-genome/JA/maker/src/../bin/chado2gff3
Installing /data/putnamlab/jillashey/pdam-genome/JA/maker/src/../bin/maker_map_ids
Installing /data/putnamlab/jillashey/pdam-genome/JA/maker/src/../bin/map_data_ids
Installing /data/putnamlab/jillashey/pdam-genome/JA/maker/src/../bin/split_gff3
Installing /data/putnamlab/jillashey/pdam-genome/JA/maker/src/../bin/AED_cdf_generator.pl
Installing /data/putnamlab/jillashey/pdam-genome/JA/maker/src/../bin/fasta_tool
Installing /data/putnamlab/jillashey/pdam-genome/JA/maker/src/../bin/genemark_gtf2gff3
Installing /data/putnamlab/jillashey/pdam-genome/JA/maker/src/../bin/ipr_update_gff
Installing /data/putnamlab/jillashey/pdam-genome/JA/maker/src/../bin/maker
Installing /data/putnamlab/jillashey/pdam-genome/JA/maker/src/../bin/maker2chado
Installing /data/putnamlab/jillashey/pdam-genome/JA/maker/src/../bin/maker2eval_gtf
Installing /data/putnamlab/jillashey/pdam-genome/JA/maker/src/../bin/maker2wap
Installing /data/putnamlab/jillashey/pdam-genome/JA/maker/src/../bin/maker2zff
Installing /data/putnamlab/jillashey/pdam-genome/JA/maker/src/../bin/map_gff_ids
Installing /data/putnamlab/jillashey/pdam-genome/JA/maker/src/../bin/tophat2gff3
Installing /data/putnamlab/jillashey/pdam-genome/JA/maker/src/../bin/train_augustus.pl
Installing /data/putnamlab/jillashey/pdam-genome/JA/maker/src/../bin/cegma2zff
Installing /data/putnamlab/jillashey/pdam-genome/JA/maker/src/../bin/fasta_merge
Installing /data/putnamlab/jillashey/pdam-genome/JA/maker/src/../bin/gff3_merge
Installing /data/putnamlab/jillashey/pdam-genome/JA/maker/src/../bin/map2assembly
Installing /data/putnamlab/jillashey/pdam-genome/JA/maker/src/../bin/zff2genbank.pl

```
Still gives me same slurm error message when I try to run, may need to: 

```
Can't locate MAKER/ConfigData.pm in @INC (you may need to install the MAKER::ConfigData module) (@INC contains: /var/spool/slurmd/job1689783/../perl/lib /var/spool/slurmd/job1689783/../lib /var/spool/slurmd/job1689783/../src/inc/perl/lib /opt/software/BioPerl/1.7.2-intel-2019b-Perl-5.30.0/lib/perl5/site_perl/5.30.0/x86_64-linux-thread-multi /opt/software/BioPerl/1.7.2-intel-2019b-Perl-5.30.0/lib/perl5/site_perl/5.30.0 /opt/software/BioPerl/1.7.2-intel-2019b-Perl-5.30.0 /opt/software/XML-LibXML/2.0132-GCCcore-8.3.0-Perl-5.30.0/lib/perl5/site_perl/5.30.0/x86_64-linux-thread-multi /opt/software/XML-LibXML/2.0132-GCCcore-8.3.0-Perl-5.30.0/lib/perl5/site_perl/5.30.0 /opt/software/XML-LibXML/2.0132-GCCcore-8.3.0-Perl-5.30.0 /opt/slurm/lib64/perl5/ /opt/software/Perl/5.30.0-GCCcore-8.3.0/lib/perl5/site_perl/5.30.0/x86_64-linux-thread-multi /opt/software/Perl/5.30.0-GCCcore-8.3.0/lib/perl5/site_perl/5.30.0 /opt/software/Perl/5.30.0-GCCcore-8.3.0/lib/perl5/5.30.0/x86_64-linux-thread-multi /opt/software/Perl/5.30.0-GCCcore-8.3.0/lib/perl5/5.30.0) at /var/spool/slurmd/job1689783/slurm_script line 41.
```

May have because my paths to files were not correct also?? Not correct above.

But very confused overall. Going to just try to do it from the basics of the basics. Using Campbell et al. 2015 paper that explained various ways to use maker. Once I do that, I'll work more into Ross code.

~~~~~~~~

#### 1) generate control files 
```
maker -CTL
```
#### 2) Edit maker opts file to specify genome, transcript, and protein path

```
# in maker_opts.ctl: 
genome=/data/putnamlab/jillashey/genome/Pdam/ReefGenomics/pdam_scaffolds.fasta
est=/data/putnamlab/jillashey/genome/Pdam/ReefGenomics/pdam_transcripts.fasta
protein=/data/putnamlab/jillashey/genome/Pdam/ReefGenomics/pdam_proteins.fasta

# remove / before beginning of path
```
#### 3) Run maker 
```
maker 2> maker.error
```

Didnt work...Gave me this error
```
STATUS: Parsing control files...
STATUS: Processing and indexing input FASTA files...
ERROR: SplitDB not created correctly

 at /opt/software/maker/3.01.03/bin/../lib/GI.pm line 1178.
        GI::split_db("/data/putnamlab/jillashey/genome/Pdam/ReefGenomics/pdam_prote"..., "protein", 10, "/data/putnamlab/jillashey/pdam-genome/JA/maker_test/pdam_scaf"..., "C") called at /opt/software/maker/3.01.03/bin/maker line 531
--> rank=NA, hostname=bluewaves.uri.edu
```

SplitDB not created correctly ? 
Maybe its something to do with the reef genomics files??

```
nano maker_round0.sh

#!/bin/bash
#SBATCH --job-name="maker"
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=250GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="maker0_out_error"
#SBATCH --output="maker0_out"

module load maker/3.01.03

maker maker_exe.ctl maker_bopts.ctl maker_opts.ctl 2> maker.error

sbatch maker_round0.sh
```

No clue...gives me nonzero exit status on bluewaves. Works when I run maker on the command line without submitting a job, but that would take forever and I'd be using all the ram on the 'home' node. 

~~

Going to attempt to use repeatMasker 

```
module load RepeatMasker/4.0.9-p2-gompi-2019b-HMMER

sbatch RepeatMasker -species "Pocillopora damicornis" -dir . -a /data/putnamlab/jillashey/genome/Pdam/ReefGenomics/pdam_scaffolds.fasta

Can't locate EMBL.pm in @INC (you may need to install the EMBL module) (@INC contains: /var/spool/slurmd/job1692505 /opt/slurm/lib64/perl5/ /opt/software/Perl/5.30.0-GCCcore-8.3.0/lib/perl5/site_perl/5.30.0/x86_64-linux-thread-multi /opt/software/Perl/5.30.0-GCCcore-8.3.0/lib/perl5/site_perl/5.30.0 /opt/software/Perl/5.30.0-GCCcore-8.3.0/lib/perl5/5.30.0/x86_64-linux-thread-multi /opt/software/Perl/5.30.0-GCCcore-8.3.0/lib/perl5/5.30.0) at /var/spool/slurmd/job1692505/slurm_script line 307.


```

Putting in -e argument

```
sbatch RepeatMasker -e ncbi -species "Pocillopora damicornis" -dir . -a /data/putnamlab/jillashey/genome/Pdam/ReefGenomics/pdam_scaffolds.fasta

# same error as above 
Can't locate EMBL.pm in @INC (you may need to install the EMBL module) (@INC contains: /var/spool/slurmd/job1692506 /opt/slurm/lib64/perl5/ /opt/software/Perl/5.30.0-GCCcore-8.3.0/lib/perl5/site_perl/5.30.0/x86_64-linux-thread-multi /opt/software/Perl/5.30.0-GCCcore-8.3.0/lib/perl5/site_perl/5.30.0 /opt/software/Perl/5.30.0-GCCcore-8.3.0/lib/perl5/5.30.0/x86_64-linux-thread-multi /opt/software/Perl/5.30.0-GCCcore-8.3.0/lib/perl5/5.30.0) at /var/spool/slurmd/job1692506/slurm_script line 307.


```

Adding path to repeatMasker in software 

sbatch /opt/software/RepeatMasker/4.0.9-p2-gompi-2019b-HMMER/RepeatMasker -e ncbi -species "Pocillopora damicornis" -dir . -a /data/putnamlab/jillashey/genome/Pdam/ReefGenomics/pdam_scaffolds.fasta
