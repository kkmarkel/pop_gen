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

## Other 
### Index fasta file with samtools
Download fasta file. Check that it is bgzipped correctly. Run
```bash
samtools-faidx $FILE
```
