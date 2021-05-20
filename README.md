# Sediment Stress

Four coral species from Hawai'i and three coral species from Florida were exposed _in situ_ to varying levels of sediment load.


### Bioinformatic pipeline

![](https://raw.githubusercontent.com/JillAshey/SedimentStress/master/Images/BioinformaticPipeline.jpeg?token=APHKO3ZVWTZNRAOJLSGSRATAU3WFC)

##### Table of Contents

- Bioinf
	- HISAT2
		- HISAT2_species_genome.md - HISAT2 script for each species
		- stringtie_HISAT2_xx.md - stringTie script for each species (following hisat2)
	- STAR
		- STAR_species_genome.md - STAR script for each species
		- stringtie_STAR_xx.md - stringTie script for each species (following star)
	- BLAST.md - BLASTp for all coral proteins
	- GenomeAnnotation_Cunning_makerPipeline.md - Annotating genes with R. Cunning MAKER pipeline
	- RNASeq_pipeline_FL.md - full bioinformatics pipeline for FL corals
	- wget_genomes.md - locations to download genome and annotation files for all coral species
- Data
	- gene_GOterms_pdam_Connelly.csv - GO terms associated with Reef Genomics gene names (obtained from Connelly [github](https://github.com/michaeltconnelly/EAPSI_Pocillopora_LPS/blob/master/data/pdam_genome_genesGO.txt))
	- gene_function_pdam_Connelly.csv - functional info associated with Reef Genomics gene names (obtained from Connelly [github](https://github.com/michaeltconnelly/EAPSI_Pocillopora_LPS/blob/master/data/pdam_genome_IDInfo.gff))
	- gene_result_pdam_NCBI.csv - NCBI gene names associated with functional info and Reef Genomics gene names
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