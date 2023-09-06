# Sediment Stress

This study evaluated the transcriptomic responses of corals with different morphological characteristics in response to different types of sediment stress. Six coral species (_Acropora cervicornis_, _Montestraea cavernosa_, _Montipora capitata_, _Orbicella faveolata_, _Pocillopora acuta_, and _Porites lobata_) were assessed at two locations (Florida and Hawai’i). Floridian corals (_A. cervicornis_, _M. cavernosa_ and _O. faveolata_) were exposed to sterilized white carbonate sediment for 18 days, whereas Hawaiian corals (_M. capitata_, _P. acuta_ and _P. lobata_) were exposed to unsterilized terrigenous red soil for up to 7 days. 

### Experimental design

![](https://raw.githubusercontent.com/JillAshey/SedimentStress/master/Images/ExperimentalDesign.png)

### Bioinformatic pipeline

![](https://raw.githubusercontent.com/JillAshey/SedimentStress/master/Images/BioinformaticPipeline.png)

## Repository overview 

- [`Bioinf`](https://github.com/JillAshey/SedimentStress/tree/master/Bioinf) - contains scripts for bioinformatic analyses 
	- [`STAR`](https://github.com/JillAshey/SedimentStress/tree/master/Bioinf/STAR) - scripts for read alignment using the STAR alignment software for the Hawai’i coral species
	- [`OrthoFinder.md`](https://github.com/JillAshey/SedimentStress/blob/master/Bioinf/OrthoFinder.md) - scripts for running OrthoFinder to identify orthologous genes between species 
	- [`RNASeq_pipeline_FL.md`](https://github.com/JillAshey/SedimentStress/blob/master/Bioinf/RNASeq_pipeline_FL.md) - scripts for RNAseq analysis for Florida coral species
	- [`RNASeq_pipeline_HI.md`](https://github.com/JillAshey/SedimentStress/blob/master/Bioinf/RNASeq_pipeline_HI.md) - scripts for RNAseq analysis for Hawai’i coral species
	- [`TranscriptLengths.md`](https://github.com/JillAshey/SedimentStress/blob/master/Bioinf/TranscriptLengths.md) - script to calculate transcript lengths of each gene and protein. This information was used in the GOSeq scripts. 
	- [`wgetgenomes.md`](https://github.com/JillAshey/SedimentStress/blob/master/Bioinf/wget_genomes.md) - download information for all reference genomes used
- [`Data`](https://github.com/JillAshey/SedimentStress/tree/master/Data)
	- [`NCBI_upload`](https://github.com/JillAshey/SedimentStress/tree/master/Data/NCBI_upload) - information on NCBI sequences upload
	- [`noaa_buoy`](https://github.com/JillAshey/SedimentStress/tree/master/Data/noaa_buoy) - water temperature data from NOAA in Florida and Hawai’i during the times of the experiments
	- [`FL_sediment_metadata_raw.csv `](https://github.com/JillAshey/SedimentStress/blob/master/Data/FL_sediment_metadata_raw.csv) - metadata from the Florida experiment
	- [`HI_sediment_metadata_raw.csv `](https://github.com/JillAshey/SedimentStress/blob/master/Data/HI_sediment_metadata_raw.csv) - metadata from the Hawai’i experiment
	- [`acerv_metadata_filtered.csv `](https://github.com/JillAshey/SedimentStress/blob/master/Data/acerv_metadata_filtered.csv) - metadata from _A. cervicornis_ from the Florida experiment
 	- [`mcap_metadata_filtered.csv `](https://github.com/JillAshey/SedimentStress/blob/master/Data/mcap_metadata_filtered.csv) - metadata from _M. capitata_ from the Hawai’i experiment
	- [`mcav_metadata_filtered.csv `](https://github.com/JillAshey/SedimentStress/blob/master/Data/mcav_metadata_filtered.csv) - metadata from _M. cavernosa_ from the Florida experiment
	- [`ofav_metadata_filtered.csv `](https://github.com/JillAshey/SedimentStress/blob/master/Data/ofav_metadata_filtered.csv) - metadata from _O. faveolata_ from the Florida experiment
 	- [`pacuta_metadata_filtered.csv `](https://github.com/JillAshey/SedimentStress/blob/master/Data/pacuta_metadata_filtered.csv) - metadata from _P. acuta_ from the Hawai’i experiment
 	- [`plob_metadata_filtered.csv `](https://github.com/JillAshey/SedimentStress/blob/master/Data/plob_metadata_filtered.csv) - metadata from _P. lobata_ from the Hawai’i experiment
- [`Images`](https://github.com/JillAshey/SedimentStress/tree/master/Images) - pipeline, QC and workflow images 
- [`Output`](https://github.com/JillAshey/SedimentStress/tree/master/Output) - output from gene expression analyses. 
	- [`DESeq2`](https://github.com/JillAshey/SedimentStress/tree/master/Output/DESeq2) - differential gene expression output
		- [`acerv`](https://github.com/JillAshey/SedimentStress/tree/master/Output/DESeq2/acerv) - differential gene expression output for _A. cervicornis_ from the Florida experiment
		- [`mcap`](https://github.com/JillAshey/SedimentStress/tree/master/Output/DESeq2/mcap) - differential gene expression output for _M. capitata_ from the Hawai’i experiment
		- [`mcav`](https://github.com/JillAshey/SedimentStress/tree/master/Output/DESeq2/mcav) - differential gene expression output for _M. cavernosa_ from the Florida experiment
		- [`ofav`](https://github.com/JillAshey/SedimentStress/tree/master/Output/DESeq2/ofav) - differential gene expression output for _O. faveolata_ from the Florida experiment
		- [`pacuta`](https://github.com/JillAshey/SedimentStress/tree/master/Output/DESeq2/pacuta) - differential gene expression output for _P. acuta_ from the Hawai’i experiment
		- [`plob`](https://github.com/JillAshey/SedimentStress/tree/master/Output/DESeq2/plob) - differential gene expression output for _P. lobata_ from the Hawai’i experiment
	- [`Figs`](https://github.com/JillAshey/SedimentStress/tree/master/Output/Figs) - figure output
		- [`GOslim`](https://github.com/JillAshey/SedimentStress/tree/master/Output/Figs/GOslim) - figures for GO slim 
			- [`BP`](https://github.com/JillAshey/SedimentStress/tree/master/Output/Figs/GOslim/BP) - GO slim figures for Biological Processes
			- [`CC`](https://github.com/JillAshey/SedimentStress/tree/master/Output/Figs/GOslim/CC) - GO slim figures for Cellular Components
			- [`MF`](https://github.com/JillAshey/SedimentStress/tree/master/Output/Figs/GOslim/MF) - GO slim figures for Molecular Functions
		- [`PCA`](https://github.com/JillAshey/SedimentStress/tree/master/Output/Figs/PCA) - PCA figures of differentiall expressed genes for Florida and Hawai’i corals 
		- [`Temperature`](https://github.com/JillAshey/SedimentStress/tree/master/Output/Figs/Temperature) - temperature figures from NOAA data from Florida and Hawai’i experimental periods 
		- [`Venn`](https://github.com/JillAshey/SedimentStress/tree/master/Output/Figs/Venn) - Venn diagram figures of shared GO terms and orthogroups between species 
		- [`acerv`](https://github.com/JillAshey/SedimentStress/tree/master/Output/Figs/acerv) - PCAs, heatmaps and GO plots for all genes and DEGs for _A. cervicornis_ from the Florida experiment
		- [`mcap`](https://github.com/JillAshey/SedimentStress/tree/master/Output/Figs/mcap) - PCAs, heatmaps and GO plots for all genes and DEGs for _M. capitata_ from the Hawai’i experiment
		- [`mcav`](https://github.com/JillAshey/SedimentStress/tree/master/Output/Figs/mcav) - PCAs, heatmaps and GO plots for all genes and DEGs for _M. cavernosa_ from the Florida experiment
		- [`ofav`](https://github.com/JillAshey/SedimentStress/tree/master/Output/Figs/ofav) - PCAs, heatmaps and GO plots for all genes and DEGs for _O. faveolata_ from the Florida experiment
		- [`pacuta`](https://github.com/JillAshey/SedimentStress/tree/master/Output/Figs/pacuta) - PCAs, heatmaps and GO plots for all genes and DEGs for _P. acuta_ from the Hawai’i experiment
		- [`plob`](https://github.com/JillAshey/SedimentStress/tree/master/Output/Figs/plob) - PCAs, heatmaps and GO plots for all genes and DEGs for _P. lobata_ from the Hawai’i experiment
		- [`upset`](https://github.com/JillAshey/SedimentStress/tree/master/Output/Figs/upset) - upset plots of shared GO terms and orthogroups between species 
	- [`GOSeq`](https://github.com/JillAshey/SedimentStress/tree/master/Output/GOSeq) - gene ontology output
		- [`acerv`](https://github.com/JillAshey/SedimentStress/tree/master/Output/GOSeq/acerv) - gene ontology output for _A. cervicornis_ from the Florida experiment
		- [`mcap`](https://github.com/JillAshey/SedimentStress/tree/master/Output/GOSeq/mcap) - gene ontology output for _M. capitata_ from the Hawai’i experiment
		- [`mcav`](https://github.com/JillAshey/SedimentStress/tree/master/Output/GOSeq/mcav) - gene ontology output for _M. cavernosa_ from the Florida experiment
		- [`ofav`](https://github.com/JillAshey/SedimentStress/tree/master/Output/GOSeq/ofav) - gene ontology output for _O. faveolata_ from the Florida experiment
		- [`pacuta`](https://github.com/JillAshey/SedimentStress/tree/master/Output/GOSeq/pacuta) - gene ontology output for _P. acuta_ from the Hawai’i experiment
		- [`plob`](https://github.com/JillAshey/SedimentStress/tree/master/Output/GOSeq/plob) - gene ontology output for _P. lobata_ from the Hawai’i experiment
		- [`all.go.BP.slim_20220417.csv`](https://github.com/JillAshey/SedimentStress/blob/master/Output/GOSeq/all.go.BP.slim_20220417.csv) - gene ontology output for Biological Processes for all species 
		- [`all.go.CC.slim_20220417.csv`](https://github.com/JillAshey/SedimentStress/blob/master/Output/GOSeq/all.go.CC.slim_20220417.csv) - gene ontology output for Cellular Components for all species 
		- [`all.go.MF.slim_20220417.csv`](https://github.com/JillAshey/SedimentStress/blob/master/Output/GOSeq/all.go.MF.slim_20220417.csv) - gene ontology output for Molecular Functions for all species 
		- [`all.go.slim_20220417.csv`](https://github.com/JillAshey/SedimentStress/blob/master/Output/GOSeq/all.go.slim_20220417.csv) - gene ontology output for all species 
	- [`OrthoFinder`](https://github.com/JillAshey/SedimentStress/tree/master/Output/OrthoFinder) - OrthoFinder output
		- [`acerv`](https://github.com/JillAshey/SedimentStress/tree/master/Output/OrthoFinder/acerv) - orthofinder and DEG output for _A. cervicornis_ from the Florida experiment
		- [`mcap`](https://github.com/JillAshey/SedimentStress/tree/master/Output/OrthoFinder/mcap) - orthofinder and DEG output for _M. capitata_ from the Hawai’i experiment
		- [`mcav`](https://github.com/JillAshey/SedimentStress/tree/master/Output/OrthoFinder/mcav) - orthofinder and DEG output for _M. cavernosa_ from the Florida experiment
		- [`ofav`](https://github.com/JillAshey/SedimentStress/tree/master/Output/OrthoFinder/ofav) - orthofinder and DEG output for _O. faveolata_ from the Florida experiment
		- [`pacuta`](https://github.com/JillAshey/SedimentStress/tree/master/Output/OrthoFinder/pacuta) - orthofinder and DEG output for _P. acuta_ from the Hawai’i experiment
		- [`plob`](https://github.com/JillAshey/SedimentStress/tree/master/Output/OrthoFinder/plob) - orthofinder and DEG output for _P. lobata_ from the Hawai’i experiment
		- [`species_compare`](https://github.com/JillAshey/SedimentStress/tree/master/Output/OrthoFinder/species_compare) - files comparing each species against one another 
		- [`AllSpecies_DEGs_ortho_20221020.csv`](https://github.com/JillAshey/SedimentStress/blob/master/Output/OrthoFinder/AllSpecies_DEGs_ortho_20221020.csv) - DEGs and corresponding orthogroups for all species 
		- [`Orthogroups.tsv`](https://github.com/JillAshey/SedimentStress/blob/master/Output/OrthoFinder/Orthogroups.tsv) - raw data from orthofinder comparison for all species 
	- [`QC`](https://github.com/JillAshey/SedimentStress/tree/master/Output/QC) - multiQC output for raw and trimmed reads 
	- [`SuppTables`](https://github.com/JillAshey/SedimentStress/tree/master/Output/SuppTables) - supplementary tables for the associated manuscript 
		- [`AlignmentOutput.csv`](https://github.com/JillAshey/SedimentStress/blob/master/Output/SuppTables/AlignmentOutput.csv) - raw and clean read counts, % of reads mapped to genome 
		- [`AlignmentOutput_Averages.csv`](https://github.com/JillAshey/SedimentStress/blob/master/Output/SuppTables/AlignmentOutput_Averages.csv) - averages of raw and clean read counts, % of reads mapped to genome 
		- [`DEG_GO_Orthogroup_SummaryTable.csv`](https://github.com/JillAshey/SedimentStress/blob/master/Output/SuppTables/DEG_GO_Orthogroup_SummaryTable.csv) - summary counts of DEG, GO and orthogroups for each species 
		- [`Orthogroups_DEGs.csv`](https://github.com/JillAshey/SedimentStress/blob/master/Output/SuppTables/Orthogroups_DEGs.csv) - summary table with the number of orthogroups with DEGs and number of DEGs in orthogroups
		- [`Shared_GO_terms.csv`](https://github.com/JillAshey/SedimentStress/blob/master/Output/SuppTables/Shared_GO_terms.csv) - GO terms shared between species 
		- [`Shared_GO.xlsx`](https://github.com/JillAshey/SedimentStress/blob/master/Output/SuppTables/Shared_GO.xlsx) - GO terms shared between species in manuscript ready table 
- [`RAnalysis`](https://github.com/JillAshey/SedimentStress/tree/master/RAnalysis) - R scripts for data analyses 
	- [`acerv`](https://github.com/JillAshey/SedimentStress/tree/master/RAnalysis/acerv) - DESeq2 and GOSeq scripts for _A. cervicornis_ from the Florida experiment
	- [`gff_wrangling`](https://github.com/JillAshey/SedimentStress/tree/master/RAnalysis/gff_wrangling) - scripts for manipulating GFFs for use with STAR software
	- [`mcap`](https://github.com/JillAshey/SedimentStress/tree/master/RAnalysis/mcap) - DESeq2 and GOSeq scripts for _M. capitata_ from the Hawai’i experiment
	- [`mcav`](https://github.com/JillAshey/SedimentStress/tree/master/RAnalysis/mcav) - DESeq2 and GOSeq scripts for _M. cavernosa_ from the Florida experiment
	- [`ofav`](https://github.com/JillAshey/SedimentStress/tree/master/RAnalysis/ofav) - DESeq2 and GOSeq scripts for _O. faveolata_ from the Florida experiment
	- [`pacuta`](https://github.com/JillAshey/SedimentStress/tree/master/RAnalysis/pacuta) - DESeq2 and GOSeq scripts for _P. acuta_ from the Hawai’i experiment
	- [`plob`](https://github.com/JillAshey/SedimentStress/tree/master/RAnalysis/plob) - DESeq2 and GOSeq scripts for _P. lobata_ from the Hawai’i experiment
	- [`DEG_GO_summary.Rmd`](https://github.com/JillAshey/SedimentStress/blob/master/RAnalysis/DEG_GO_summary.Rmd) -  script summarizing DEG and GO information 
	- [`GOslim.Rmd`](https://github.com/JillAshey/SedimentStress/blob/master/RAnalysis/GOslim.Rmd) - script assigning GO terms into GO slim bins 
	- [`NOAA_temperatures`](https://github.com/JillAshey/SedimentStress/blob/master/RAnalysis/NOAA_temperatures.Rmd) - script plotting NOAA temperatures from Florida and Hawai’i
	- [`Orthogroups.Rmd`](https://github.com/JillAshey/SedimentStress/blob/master/RAnalysis/Orthogroups.Rmd) - script analyzing orthogroups in the context of DEGs for each species 
	- [`PCA_plots.Rmd`](https://github.com/JillAshey/SedimentStress/blob/master/RAnalysis/PCA_plots.Rmd) - script to make PCA figure for manuscript 
	- [`RNASeqPower.Rmd`](https://github.com/JillAshey/SedimentStress/blob/master/RAnalysis/RNASeqPower.Rmd) - script to run power analysis on RNASeq data





