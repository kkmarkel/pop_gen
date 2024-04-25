# ADMIXTURE
## About
ADMIXTURE is a program for estimating ancestry in a model-based manner from large autosomal SNP genotype datasets, where the individuals are unrelated (for example, the individuals in a case-control association study).  
ADMIXTURE’s input is binary PLINK (.bed), ordinary PLINK (.ped), or EIGENSTRAT (.geno) formatted files and its output is simple space-delimited files containing the parameter estimates.  
To use ADMIXTURE, you need an input file and an idea of K, your belief of the number of ancestral populations. You should also have the associated support files alongside your main input file, in the same directory. For example, if your primary input file is a .bed file, you should have the associated .bim (binary marker information file) and .fam (pedigree stub file) files in the same directory. If your primary input file is a .ped or .geno file, a corresponding PLINK style .map file should be in the same directory.  
## Installation
```bash
conda create -n admixture
conda install mamba
mamba install bioconda::plink
mamba install bioconda::admixture
```
## Usage

### VCF to bed-bim-fam
```bash
# Generate the input file in plink format
plink --vcf /home/data/vcf/$FILE.vcf.gz --make-bed --out $FILE --allow-extra-chr
```
### LD pruning
```bash
plink --bfile $FILE --indep-pairwise 50 10 0.1
# (output indicating number of SNPs targeted for inclusion/exclusion)

plink --bfile $FILE --extract plink.prune.in --make-bed --out prunedData
```
### Run ADMIXTURE
>[!note]- 
> ADMIXTURE can run with .bed, .ped, or .geno files
```bash
admixture prunedData.bed 6
```
ADMIXTURE starts running. When it finishes, it will output some estimates:
```bash
$ ls
prunedData.bed prunedData.bim prunedData.fam prunedData.6.Q prunedData.6.P
```
There is an output file for each parameter set: **Q (the ancestry fractions)**, and **P (the allele frequencies of the inferred ancestral populations)**  
**K**-- number of populations that was assumed for the analysis  
You can run ADMIXTURE on a input file located in another directory, but its output will
always be put in the current working directory.  

If you also wanted standard errors:
```bash
admixture -B prunedData.bed 6
```
This will perform point estimation and then will also use a bootstrapping procedure to
calculate the standard errors. Note that (point-estimation & bootstrapping) takes considerably longer than point-estimation alone, so you will have to be patient. Eventually it
will finish, yielding point estimates and standard errors:
```bash
% ls
prunedData.6.Q prunedData.6.P prunedData.6.Q_se
```
The “_se” file is in the same unadorned file format as the point estimates.

If your analyses are taking a long time, first consider why.   
* Is your dataset huge? Do
you really need to analyze all the markers, or can you thin the marker set ([LD_pruning?](#ld-pruning))
* Larger K => larger time
* Consider using multithreading
```bash
admixture prunedData.bed 6 -j4
# for four cpu
```
## Links
[Software manual 1.3.0](https://dalexander.github.io/admixture/admixture-manual.pdf)