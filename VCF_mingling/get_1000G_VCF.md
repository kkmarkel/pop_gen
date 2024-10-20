# Download and process 1000G data to get a VCF file

## Download
Information about 1000G data: https://www.internationalgenome.org/data-portal/data-collection

To download `1000 Genomes 30x on GRCh38` data
```bash
for chr in {1..22}; do \
    wget "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"  "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi" ; 
done
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz 
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz.tbi
```

## Process and convert the 1000 Genomes files to BCF
* Ensure that multi-allelic calls are split and that indels are left-aligned compared to reference genome (1st pipe)
* Set the ID field to a unique value: CHROM:POS:REF:ALT (2nd pipe)
* Remove duplicates (3rd pipe)

`-I +'%CHROM:%POS:%REF:%ALT'` means that unset IDs will be set to `CHROM:POS:REF:ALT`

`-x ID -I +'%CHROM:%POS:%REF:%ALT'` first erases the current ID and then sets it to `CHROM:POS:REF:ALT`

```bash
# very long!
for chr in {1..22}; do
    bcftools norm -m-any --check-ref w -f /mnt/d/ref/gencode/GRCh38.primary_assembly.genome.fa.gz \
    1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz | \
    bcftools annotate -x ID -I +'%CHROM:%POS:%REF:%ALT' | \
    bcftools norm --rm-dup both \
    -Ob -o 1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.bcf ;

    bcftools 1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.bcf ;
done
```