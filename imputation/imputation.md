# Imputation workflow
## References
* https://www.protocols.io/run/genotype-imputation-workflow-v3-0-xbgfijw?
## Requirements and preparatory steps
### Conda environment 
```bash
mamba env create -f imputation.yml
conda activate imputation
```
### Reference genome and genetic map files
### Imputation reference panel files
For increased imputation accuracy, it's recommended to use a population-specific imputation reference panel, if available.  
If population-specific reference data is not available, for instance 1000 Genomes Project (1000 GP) (www.nature.com/articles/nature15393) data can be used instead.  

GRCh37/hg19 files are available at Beagle site already processed to be compatible with Beagle:
http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.bref/ 
#### 1000 genomes phase 3 High Coverage SNP+Indel
```bash
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr{{1..22},X}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz{,.tbi}

# for chrX the name is a little different, so:
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz{,.tbi}
# after that rename it accordingly
```
##### Minimum QC
For 1000GP data:

1. Remove the rare variants, here singletons and doubletons by setting AC threshold with 'bcftools view'.
2. Split multiallelic sites to biallelic records with 'bcftools norm'.
3. Keep only SNPs and INDELs with 'bcftools view'. Here, the 1000GP data included a tag VT in the INFO field and data contain also structural variants which should be excluded.
5. Align the variants to reference genome with 'bcftools norm' in order to have the REF and ALT alleles in the shortest possible representation and to confirm that the REF allele matches the reference genome, additionally remove duplicate variants (-d none).
6. After alignment, remove multiallelic records with 'bcftools view', since these are formed during the alignment if the REF does not match with the reference genome.
7. Finally, remove sites containing missing data with 'bcftools view'.

```bash
# Generate a chromosome renaming file
for CHR in {1..23} X ; do 
    echo ${CHR} chr${CHR}
done >> chr_names.txt

# min_qc.sh script
# Multiple processing commands piped together
for CHR in {1..22} X; do
    bcftools annotate --rename-chrs chr_names.txt \
        ALL.chr${CHR}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz -Ou | \
    bcftools view -e 'INFO/AC<3 | INFO/AN-INFO/AC<3' -Ou | \
    bcftools norm -m -any -Ou | \
    bcftools view -i 'INFO/VT="SNP" | INFO/VT="INDEL"' -Ou | \
    bcftools norm -f /storage/kkmarkel/pop/ref/gencode/GRCh38.primary_assembly.genome.fa.gz -d none -Ou | \
    bcftools view -m 2 -M 2 -Ou | \
    bcftools view -g ^miss -Oz -o 1000GP_chr${CHR}.vcf.gz
done
```
After this there is a possible step to preserve multiallelic sites - setting unique IDs.
```bash
for CHR in {1..23}; do
    bcftools annotate \
    --set-id '%CHROM\_%POS\_%REF\_%ALT' \
    1000GP_chr${CHR}.vcf.gz  \
    -Oz -o 1000GP_SNPID_chr${CHR}.vcf.gz
done
```


#### 1000 genomes 20190312 Biallelic SNP+INDEL (N=2548)
```bash
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr{{1..22},X}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz{,.tbi}
```
> ![note]- Important note
> The reference panel files should contain: 
> * phased genotypes,
> * chromosome names with 'chr' and chromosome X as 'chrX', 
> * all variants as biallelic records,
> * only SNPs and INDELs,
> * only unique variants, 
> * non-missing data, 
> * chrX as diploid genotypes, and
> * only unique IDs

##### Minimum QC
For 1000GP data:
1. Generate a text file containing space-spearated old and new chromosome names. This is required to rename the numerical chromosome names with 'chr' tag. Apply the new chromosome names with 'bcftools annotate'.
2. Remove the rare variants, here singletons and doubletons by setting AC threshold with 'bcftools view'.
3. Split multiallelic sites to biallelic records with 'bcftools norm'.
4. Keep only SNPs and INDELs with 'bcftools view'. Here, the 1000GP data included a tag VT in the INFO field and data contain also structural variants which should be excluded.
5. Align the variants to reference genome with 'bcftools norm' in order to have the REF and ALT alleles in the shortest possible representation and to confirm that the REF allele matches the reference genome, additionally remove duplicate variants (-d none).
6. After alignment, remove multiallelic records with 'bcftools view', since these are formed during the alignment if the REF does not match with the reference genome.
7. Finally, remove sites containing missing data with 'bcftools view'.

```bash
# Generate a chromosome renaming file
for CHR in {1..23} X ; do 
    echo ${CHR} chr${CHR}
done >> chr_names.txt

# min_qc.sh script
# Multiple processing commands piped together
for CHR in {1..22} X; do
    bcftools annotate --rename-chrs chr_names.txt \
        ALL.chr${CHR}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz -Ou | \
    bcftools view -e 'INFO/AC<3 | INFO/AN-INFO/AC<3' -Ou | \
    bcftools norm -m -any -Ou | \
    bcftools view -i 'INFO/VT="SNP" | INFO/VT="INDEL"' -Ou | \
    bcftools norm -f /storage/kkmarkel/pop/ref/gencode/GRCh38.primary_assembly.genome.fa.gz -d none -Ou | \
    bcftools view -m 2 -M 2 -Ou | \
    bcftools view -g ^miss -Oz -o 1000GP_chr${CHR}.vcf.gz
done
```
After this there is a possible step to preserve multiallelic sites - setting unique IDs.
```bash
for CHR in {1..23}; do
    bcftools annotate \
    --set-id '%CHROM\_%POS\_%REF\_%ALT' \
    1000GP_chr${CHR}.vcf.gz  \
    -Oz -o 1000GP_SNPID_chr${CHR}.vcf.gz
done
```
#TODO include that into the previous script?

##### Convert haploid genotypes to homozygous diploids
Often chrX is represented as haploid genotypes for males, however, Beagle can only handle diploid genotypes.   
The command here produces unphased diploid genotypes. But since the haploid genotypes are in diploid format as REF/REF or ALT/ALT, we can simply set the phase for those alleles with a simple sed replacement. 

The chrX ploidy can be corrected as follows:
```bash
# Fix the chromosome X ploidy to phased diploid
# Requires a ploidy.txt file containing 
# space-separated CHROM,FROM,TO,SEX,PLOIDY 
echo "chrX 1 156040895 M 2" > ploidy.txt
bcftools +fixploidy \
    1000GP_SNPID_chrX.vcf.gz -Ov -- -p ploidy.txt | \
    sed 's#0/0#0\|0#g;s#1/1#1\|1#g' | \
bcftools view -Oz -o 1000GP_chr23.vcf.gz
```
##### Duplicate ID removal
Remove duplicate IDs. If you wish to preserve all multiallelic sites, replace the ID column with a unique ID e.g. CHR_POS_REF_ALT (as indicated in [Minimum QC](#minimum-qc)).  
```bash
for CHR in {1..23}; do
    bcftools query -f '%ID\n' 1000GP_SNPID_chr${CHR}.vcf.gz | \
    sort | uniq -d > 1000GP_SNPID_chr${CHR}.dup_id

    if [[ -s 1000GP_SNPID_chr${CHR}.dup_id ]]; then
    	bcftools view -e ID=@1000GP_SNPID_chr${CHR}.dup_id \
    	1000GP_SNPID_chr${CHR}.vcf.gz \
        -Oz -o 1000GP_filtered_chr${CHR}.vcf.gz
    else 
    	mv 1000GP_SNPID_chr${CHR}.vcf.gz \
        1000GP_filtered_chr${CHR}.vcf.gz
    fi
done
```
##### Reference panel allele frequencies
To be continued ...
