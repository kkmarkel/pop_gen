# VCF mingling tasks using bcftools

## Left-align and normalize 
```bash
bcftools norm file.vcf -f ref.fa -m -any --multi-overlaps 0 -Oz -o file.norm.vcf.gz
```
## Remove indels
```bash
bcftools view file.vcf --type snps -Oz -o file.snp.vcf.gz
```

## Fix alt and ref
>[!note] Reference allele mismatches are a huge pain in the ass. 
```bash
# checking the reference allele mismatches
bcftools +fixref test.bcf -- -f ref.fa
# or 
bcftools norm --check-ref e -f /path/to/reference.fasta input.vcf.gz


-m flip-all
#flipping any sites, including ambiguous A/T and C/G pairs

-m id
# flip based on dbSNP rs ids
```
#TODO
- [] Does `bcftools norm ... --check-ref s` work in the same way as `bcftools +fixref ... -m flip-all`?

## Extract the proportion of missing data per sample in a VCF/BCF file
Simply replace `<<FILE>>` with your properly formated VCF/BCF file name (2 places).
Required bcftools v. 1.2+.

```shell
paste \
<(bcftools query -f '[%SAMPLE\t]\n' <<FILE>> | head -1 | tr '\t' '\n') \
<(bcftools query -f '[%GT\t]\n' <<FILE>> | awk -v OFS="\t" '{for (i=1;i<=NF;i++) if ($i == "./.") sum[i]+=1 } END {for (i in sum) print i, sum[i] / NR }' | sort -k1,1n | cut -f 2)
```

## Liftover
Using the build from `score` repo
```bash
# docker image from this bundle https://software.broadinstitute.org/software/mochawdl/
docker pull us.gcr.io/mccarroll-mocha/bcftools:1.17-20230919
```

## Merge VCF files
Create text file with VCF file names, one file name per line
```bash
bcftools merge -l vcf.txt -m none -Ou | bcftools sort - -Oz -o file.vcf.gz
```

## Other 
### Index fasta file with samtools
Download fasta file. Check that it is bgzipped correctly. Run
```bash
samtools-faidx $FILE
```
