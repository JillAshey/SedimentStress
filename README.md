# Sediment Stress

This study evaluated the transcriptomic responses of corals with different morphological characteristics in response to different types of sediment stress. Six coral species (_Acropora cervicornis_, _Montestraea cavernosa_, _Montipora capitata_, _Orbicella faveolata_, _Pocillopora acuta_, and _Porites lobata_) were assessed at two locations (Florida and Hawai’i). Floridian corals (_A. cervicornis_, _M. cavernosa_ and _O. faveolata_) were exposed to sterilized white carbonate sediment for 18 days, whereas Hawaiian corals (_M. capitata_, _P. acuta_ and _P. lobata_) were exposed to unsterilized terrigenous red soil for up to 7 days. 

### Experimental design

![](https://raw.githubusercontent.com/JillAshey/SedimentStress/master/Images/ExperimentalDesign.png)

### Bioinformatic pipeline

![](https://raw.githubusercontent.com/JillAshey/SedimentStress/master/Images/BioinformaticPipeline.png)

## Repository overview 

- [`Bioinf`](https://github.com/JillAshey/SedimentStress/tree/master/Bioinf) - contains scripts for bioinformatic analyses 
	- [`STAR`](https://github.com/JillAshey/SedimentStress/tree/master/Bioinf/STAR) - scripts for read alignment using the STAR alignment software for the Hawai’i coral species
		- [`mcap_STAR.md`](https://github.com/JillAshey/SedimentStress/blob/master/Bioinf/STAR/mcap_STAR.md) - script for read alignment of _M. capitata_ reads
		- [`pacuta_STAR.md`](https://github.com/JillAshey/SedimentStress/blob/master/Bioinf/STAR/pacuta_STAR.md) - script for read alignment of _P. acuta_ reads
		- [`plob_STAR.md`](https://github.com/JillAshey/SedimentStress/blob/master/Bioinf/STAR/plob_STAR.md) - script for read alignment of _P. lobata_ reads
	- [`OrthoFinder.md`](https://github.com/JillAshey/SedimentStress/blob/master/Bioinf/OrthoFinder.md) - scripts for running OrthoFinder to identify orthologous genes between species 
	- [`RNASeq_pipeline_FL.md`](https://github.com/JillAshey/SedimentStress/blob/master/Bioinf/RNASeq_pipeline_FL.md) - scripts for RNAseq analysis for Florida coral species
	- [`RNASeq_pipeline_HI.md`](https://github.com/JillAshey/SedimentStress/blob/master/Bioinf/RNASeq_pipeline_HI.md) - scripts for RNAseq analysis for Hawai’i coral species
	- [`TranscriptLengths.md`](https://github.com/JillAshey/SedimentStress/blob/master/Bioinf/TranscriptLengths.md) - script to calculate transcript lengths of each gene and protein. This information was used in the GOSeq scripts. 
	- [`wgetgenomes.md`](https://github.com/JillAshey/SedimentStress/blob/master/Bioinf/wget_genomes.md) - download information for all reference genomes used
- [`Data`](https://github.com/JillAshey/SedimentStress/tree/master/Data)