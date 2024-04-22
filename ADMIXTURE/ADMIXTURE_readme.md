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
## Links
[Software manual 1.3.0](https://dalexander.github.io/admixture/admixture-manual.pdf)