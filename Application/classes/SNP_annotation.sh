<<COMMENT
Name:            
Author:          manuel.dominguezbecerra@nhs.net
Description:     
States:          RELEASE
COMMENT
#!/bin/sh 

# Clinicians request delete variants that don't pass 
# all filters. So, get only PASS variants.
# we can do this with BCFtools
module load bio/BCFtools/1.9-foss-2019b
bcftools view -f PASS /re_gecip/machine_learning/LOEUF_classifier_tool/Data_output/tmp/filtered_SNPs_1.vcf > /re_gecip/machine_learning/LOEUF_classifier_tool/Data_output/tmp/filtered_SNPs_2.vcf

# This could be also done by grep (Future DEV)

# Gel documentation says that VEP module doesn't work
# if there is another module loaded when using VEP module
# So, we need to unload them (I have seen that this is not true as
# as the app works without this but just in case)
# This consumes by the way half of the execution time of this app.
module purge
 
# Now we can load VEP module
module load bio/VEP/99.1-foss-2019a-Perl-5.28.1


# Run the "export" line in order to add the relevant version of LOFTEE to the Perl library variable
# (note: you need to choose GRCh37 or GRCh38 here because, the code for processing GRCh37 or GRCh38 
# data is found at different locations in Gel ) - this is required by the LOFTEE plugin:
export PERL5LIB=$PERL5LIB:${LOFTEE38}

# Run VEP + LOFTEE, and using the variables that are relevant to our case in the  "--plugin"  option line 
vep \
    --input_file /re_gecip/machine_learning/LOEUF_classifier_tool/Data_output/tmp/filtered_SNPs_2.vcf \
    --format vcf \
    --offline \
    --species homo_sapiens \
    --cache \
    --cache_version 99 \
    --assembly GRCh38 \
    --dir_cache ${CACHEDIR} \
    --verbose \
    --no_stats \
    --fasta /public_data_resources/reference/GRCh38/GRCh38Decoy_no_alt.fa \
    --plugin CADD,/public_data_resources/CADD/v1.5/GRCh38/whole_genome_SNVs.tsv.gz \
    --plugin LoF,loftee_path:${LOFTEE38},human_ancestor_fa:${LOFTEE38HA},gerp_bigwig:${LOFTEE38GERP},conservation_file:${LOFTEE38SQL} \
    --force_overwrite \
    --fields ZYG \
    --individual all \
    --gencode_basic \
    --hgvs \
    --output_file /re_gecip/machine_learning/LOEUF_classifier_tool/Data_output/tmp/filtered_and_annotated_VEP_SNPs.vcf

#  Take high confident LoF variants (Lof=HC)
grep -E 'LoF=HC' /re_gecip/machine_learning/LOEUF_classifier_tool/Data_output/tmp/filtered_and_annotated_VEP_SNPs.vcf > /re_gecip/machine_learning/LOEUF_classifier_tool/Data_output/tmp/filtered_and_annotated_VEP_HC_SNPs.vcf

# Load again modules we need
# this is in case more than 1 genome is analysed.
module load bio/BEDTools/2.19.1
module load bio/BCFtools/1.9-foss-2019b
