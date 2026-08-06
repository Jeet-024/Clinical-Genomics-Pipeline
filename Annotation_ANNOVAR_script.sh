#This script is written to execute ANNOVAR annotation tool (which is a perl based program and requires databases to be downloaded locally or in server which ever in use)
# the protocol flag mentions the databases used to annotate the variants
#!/bin/bash

set -eo pipefail

REF="/path/to/human/genome/reference/hg38.fa"
bcftools norm -m -both -f "$REF" input_vcf.gz -Oz -o output_norm_vcf.gz #This is a crucial step before annotation can be done using ANNOVAR, it splits any multiallelic sites into biallelic for non-erroneous annotation

perl /mnt/lab_data/tools/annovar/table_annovar.pl /path/to/the/output_norm_vcf.gz /mnt/lab_data/tools/annovar/humandb \
-buildver hg38 \
-outfile /path/to/output/file/that/will/be/created \
-remove \
-protocol refGene,clinvar_20240917,dbnsfp47a,gnomad41_exome,gnomad41_genome,indigen,ALL.sites.2015_08,AFR.sites.2015_08,AMR.sites.2015_08,EAS.sites.2015_08,EUR.sites.2015_08,SAS.sites.2015_08 \
-operation g,f,f,f,f,f,f,f,f,f,f,f \
-otherinfo -nastring . -vcfinput -polish

echo "Annotation done successfully"
