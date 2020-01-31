# MEL long read sequencing analysis
This repository contains code to process and analyze PacBio long read sequencing data related to Reimer et al., 2020. Raw and processed data associated with this publication can be found at GEO with accession ID [GSE144205](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144205).

## Contents
1) Code related to filtering and mapping data
2) Code related to quantifying splicing

### raw_data_processing
* [data](./raw_data_processing/data) folder contains:
	* genome file of GLOBE vector sequence used to map HBB-targeted long read sequencing data
	* genome index of GLOBE vector sequene
	
* [scripts](./raw_data_processing/scripts) folder contains:
	* python script for removing polyA reads from mm10 genome wide data (filter_polyA.py)
	* python script for removing polyA reads from HBB-targeted data (filter_polyA_HBB.py)


### Splicing Quantification
Custom scripts used to quantify splicing status are found in the [splicing_quantification](MEL_LRS/splicing_quantification) directory.

