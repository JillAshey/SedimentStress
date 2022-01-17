# Sediment Stress

This study evaluated the transcriptomic responses of corals with different morphological characteristics in response to different types of sediment stress. Five coral species (_Acropora cervicornis_, _Montestraea cavernosa_, _Orbicella faveolata_, _Pocillopora acuta_, and _Porites lobata_) were assessed at two locations (Florida and Hawai’i). Floridian corals (_A. cervicornis_, _M. cavernosa_ and _O. faveolata_) were exposed to sterilized white carbonate sediment for 30 days, whereas Hawaiian corals (_P. acuta_ and _P. lobata_) were exposed to unsterilized terrigenous red soil for up to 7 days. 

### Bioinformatic pipeline

![](https://raw.githubusercontent.com/JillAshey/SedimentStress/master/Images/BioinformaticPipeline.jpeg?token=APHKO3ZVWTZNRAOJLSGSRATAU3WFC)

##### Table of Contents

- Bioinf
	- STAR
		- STAR_species_genome.md - STAR script for each species
		- stringtie_STAR_xx.md - stringTie script for each species (following star)
	- BLAST.md - BLASTp for all coral proteins
	- RNASeq_pipeline_FL.md - full bioinformatics pipeline for FL corals
	- wget_genomes.md - locations to download genome and annotation files for all coral species
- Data
	- metadata_species_raw_filtered.csv - metadata associated with species and sample names 
	- sediment_FL_metadata.csv - FL metadata for all species
	- sediment_HT_metadata.csv - HI metadata for all species 
	- universal_pdam_Connelly.csv - list of Reef Genomics gene names 
- Images
	- QC
		- raw
			- FL - pngs of multiQC results for FL raw samples 
			- HI - pngs of multiQC results for FL raw samples 
		- trimmed 
			- FL - pngs of multiQC results for FL trimmed samples 
			- HI - pngs of multiQC results for FL trimmed samples 
- Output
	- DESeq2
	- GOSeq
	- Plots
	- QC
- RAnalysis