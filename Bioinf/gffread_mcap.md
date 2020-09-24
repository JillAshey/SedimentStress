Using gffread to convert mcap gff file to gtf file. May be more compatible with STAR

```
nano gffread_mcap.sh


#!/bin/bash
#SBATCH --job-name="gffread"
#SBATCH -t 48:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jillashey@uri.edu
#SBATCH --account=putnamlab
#SBATCH --error="gffread_mcap_out_error"
#SBATCH --output="gffread_mcap_out"

module load gffread/0.12.2-GCCcore-8.3.0 

cd /data/putnamlab/jillashey/genome/Mcap

gffread Mcap.GFFannotation.fixed.gff -T -o mcap.annotation.gtf

sbatch gffread_mcap.sh

```

Submitted batch job 1761614
Success! Really quick too

```
gffread Pcomp.GFFannotation.fixed_transcript.gff -T -o pcomp.annotation.gtf

```