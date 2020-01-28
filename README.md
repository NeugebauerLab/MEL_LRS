# MEL_LRS
This repository contains code to process and analyze PacBio long read sequencing data related to Reimer et al., 2020. Raw and processed data associated with this publication can be found at GEO with accession ID [GSE144205](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144205).

## Contents
1) Code to pre-process, map, and filter raw data
2) Code to generate splicing and readthrough quantification files
3) Code to generate figures related to long read sequencing in Reimer et al., 2020


### Pre-processing, mapping, and filtering raw data
The easiest way to download the raw FASTQ data associated with this project is from the [Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra). A snakemake workflow is provided to trim, map, and filter the data for both genome-wide and targeted datasets. Download the [raw_data_processing](/raw_data_processing) directory, and move raw data in FASTQ format to a new folder called ```0_samples```. You will also need to add to this directory a genome file (in FASTA format) and an index file for mapping to the mm10 genome. These files should be called ```genome.fa``` and ```genome.fai```. Information for downloading the mm10 genome sequence can be found at [UCSC](https://hgdownload.soe.ucsc.edu/downloads.html). 

1. Create conda environment with necessary dependencies from conda environment file:
```
conda env create -f environment.yml
```
2. Activate conda environment:
```
conda activate snakemake
```
3. Run snakemake command for processing either genome-wide (```Snakefile_mm10```) or targeted (```Snakefile_HBB```) data:
```
snakemake -s Snakefile_*
```

Note that Porechop is used to remove sequencing adapters and de-concatenate reads. Instructions for Porechop installation are [here](https://github.com/rrwick/Porechop). In order to specify the adapters that need to be trimmed, you must edit the ```adapters.py``` file so that it contains the adapters used in this dataset:
```
ADAPTERS = [
        Adapter('SMARTER-adapters',
                start_sequence=('SMARTerIIA_5prime', 'AAGCAGTGGTATCAACGCAGAGTAC'),
                end_sequence=('SMARTerIIA_3prime', 'GTACTCTGCGTTGATACCACTGCTT')),
					
		    Adapter('tLRS-adapters-fwd',
		            start_sequence=('common-seq-fwd', 'GACGTGTGCTCTTCCGATCT'),
		            end_sequence=('IIA-seq-fwd', 'GTACTCTGCGTTGATACCACTGCTT')),
					
		    Adapter('tLRS-adapters-rev',
		            start_sequence=('common-seq-rc', 'AGATCGGAAGAGCACACGTC'),
		            end_sequence=('IIA-seq-rrc', 'AAGCAGTGGTATCAACGCAGAGTAC'))
                ]
```
