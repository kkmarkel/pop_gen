# RFMix
## Installation
```bash
mamba env create -f rfmix.yml
conda activate rfmix
```
## Gene map
Gene maps premade for Eagle are available for download from https://alkesgroup.broadinstitute.org/Eagle/downloads/tables/
```bash
wget https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/tables/genetic_map_hg38_withX.txt.gz
gunzip genetic_map_hg38_withX.txt.gz
```
## Run
```bash
rfmix -f <query VCF/BCF file>
	-r <reference VCF/BCF file>
	-m <sample map file>
	-g <genetic map file>
	-o <output basename>
	--chromosome=<chromosome to analyze> 
```