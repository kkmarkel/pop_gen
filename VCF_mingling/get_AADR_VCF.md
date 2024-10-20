# Download and process AADR data to get a VCF file

## Download
Information about AADR (The Allen Ancient DNA Resource): A curated compendium of ancient human genomesdata
* [old site](https://reich.hms.harvard.edu/allen-ancient-dna-resource-aadr-downloadable-genotypes-present-day-and-ancient-dna-data)
* [new site](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/FFIDCW)

```bash
cd /mnt/d/data/AADR
```
Get 1240K+HO (HO=Human Origins) data - select the needed files on Dataverse and run download.

## Convert from Eigenstrat to PACKEDPED
Create a parfile `par.EIGENSTRAT.PACKEDPED.1240K_HO.txt` with the following contents:   
```bash
genotypename:    v62.0_HO_public.geno
snpname:         v62.0_HO_public.snp
indivname:       v62.0_HO_public.ind
outputformat:    PACKEDPED
genotypeoutname: v62.0_HO_public.bed
snpoutname:      v62.0_HO_public.bim
indivoutname:    v62.0_HO_public.fam
```
Run:
```bash
mamba activate geno
convertf -p par.EIGENSTRAT.PACKEDPED.1240K_HO.txt
```
DO NOT use just this packedped - it can have issues with SNP order and reference allele choices. I decided to correct it with bcftools later.

## Convert from PACKEDPED to VCF

```bash
plink2 \
--bfile v62.0_HO_public \
--export vcf-4.2 id-paste=iid \
--fa /mnt/d/ref/hg19.ensembl.111/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz \
--ref-from-fa  \
--geno 1.0 \
--not-chr 0,M,XY \
--out v62.0_HO_public
# 21945 samples (8385 females, 12741 males, 819 ambiguous; 21945 founders)
# 584131 variants loaded from v62.0_HO_public.bim.
bgzip v62.0_HO_public.vcf
tabix v62.0_HO_public.vcf.gz
```

## Normalize VCF and convert to BCF
* Ensure that multi-allelic calls are split and that indels are left-aligned compared to reference genome (1st pipe)
* Sets the ID field to a unique value: CHROM:POS:REF:ALT (2nd pipe)  
`-x ID -I +'%CHROM:%POS:%REF:%ALT'` first erases the current ID and then sets it to `CHROM:POS:REF:ALT`
```bash
bcftools norm v62.0_HO_public.vcf.gz -m-any --check-ref ws -f /mnt/d/ref/hg19.ensembl.111/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz | \
bcftools annotate -x ID -I +'%CHROM:%POS:%REF:%ALT' -Ob -o v62.0_HO_public.bcf


```