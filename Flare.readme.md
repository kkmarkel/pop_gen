# Flare 
## Installation
Requires Java 1.8 or later  
```bash
wget https://faculty.washington.edu/browning/flare.jar
# download the executable
```
## Usage
```bash
java -jar flare.jar 
# print a summary of the command line arguments 
```
to run Flare:
```bash
java -Xmx[GB]g -jar flare.jar [arguments]
```
where **[GB]** is the maximum number of gigabytes of memory to use, and **[arguments]** is a space-separated list of parameter values, each expressed as parameter=value.
### Required parameters
The flare program has five required parameters. Two of the required parameters specify Variant Call Format (VCF) files. A VCF record may have multiple ALT alleles and must include a genotype (GT) FORMAT subfield. **All genotypes must be phased and have no missing alleles.**  
If a VCF file has unphased or missing genotypes, you can phase the genotypes and fill in the missing genotypes using the Beagle program. Any input file with a name ending in ".gz" is assumed to be gzip-compressed.
## Links
* [Github page](https://github.com/browning-lab/flare)