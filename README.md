# MEL long read sequencing analysis
This repository contains code to process and analyze PacBio long read sequencing data related to Reimer et al., 2020. Raw and processed data associated with this publication can be found at GEO with accession ID [GSE144205](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144205).

## Contents
1) Code related to filtering and mapping data
2) Code related to quantifying splicing

### Raw Data Processing
* [data](./raw_data_processing/data) folder contains:
	* genome file of GLOBE vector sequence used to map HBB-targeted long read sequencing data
	* genome index of GLOBE vector sequene
	
* [scripts](./raw_data_processing/scripts) folder contains:
	* filter_polyA.py: script for removing polyA reads from mm10 genome wide data
	* filter_polyA_HBB.py: script for removing polyA reads from HBB-targeted data


### Splicing Quantification
* [splicing_quantification](MEL_LRS/splicing_quantification) directory folder contains:
	* calculate_distance_to_splicing.py: script to calculate the last block size and upstream intron coordinates from a BED12 file
	* classify_cigar.py: script to separate a SAM file based on whether CIGAR string indicates spliced or unspliced
	* classify_intron_splicing_status.py: script to classify splicing status of each intron-spanning region of long read sequences in BED12 format
