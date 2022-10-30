'''
Name:            short_variant.py
Author:          manuel.dominguezbecerra@nhs.net
Description:     All SNPs and indel manipulation is carried out in this script
States:          RELEASE
'''

import os
import io
import sys
import argparse
import subprocess
import pandas as pd
import numpy as np 
from classes import log
from classes import import_data
from classes.read_file import ReadFile
from classes import data_presentation2


class Snp_manipulation:
    def __init__(self,Standard_VCF_path, Structural_VCF, Genes, file_name):
        self.Standard_VCF_path = Standard_VCF_path
        self.Structural_VCF = Structural_VCF
        self.Genes = Genes
        self.file_name = file_name

    def Identify_affected_transcripts(self,Standard_VCF_path=None, Structural_VCF=None,Genes=None, file_name=None):
        '''
        Intead of filtering SNPs overlapping CNVs and genes
        They wants all SNPs that overlap with the exons of the affected tranacripts
        This is a request by the clinical team
        To do this, we need unique transcripts ID
        Then, we take their start and end coordinates
        Fortunately, I got this data from the dataset when annotation was done
        Finally, we generate the bed file we need for the filtering process

        GEL is not perfect and I have found that this not always have the counterpart standard VCF file
        in the same directory or was not created in that partciular version or the file name is different.

        I have also identified that some  standard VCF files are not well formed and bedtools return 
        "Error: Unable to open file or unable to determine types for file".
        (I think this is because the standard VCF contains all snp and indel of genome, making the file very large)

        If this occurs, I have some try/except or if/else lines to detect this and stop the app.
        CNV were annotated but genes are not annotated because this is done in a downstream function.
        I have created a alternative data_presentation2 script that annotates genes in the event of 
        not finding HC LoF variants, not finding the standard VCF file or when bedtools return this error.
        The outpout generated is identical that the one with SNPs and indels but logically columns
        that contain short variant info are empty. This worked in the 71k genome analysis carried out.   
        '''
        # First, let check that the Standar CVF file exists

        Standard_VCF_path=self.Standard_VCF_path
        genes=self.Genes
        file_name=self.file_name
        if os.path.exists(Standard_VCF_path):
                log.log_messages('INFO:        Standard VCF file exists')
        else:
                message = 'ERROR: ' +  Standard_VCF_path +  'is not found!' 
                log.log_messages(message)
                log.log_messages('INFO:         Annotating genes by alternative route')
                data_presentation2.hand_snp_annotation_error(genes,file_name)

        # Get our CNVs
        Structural_VCF=self.Structural_VCF

        # Delete CNVs that doesnt affect genes or exons
        Structural_VCF = Structural_VCF[~Structural_VCF['Constitutive_exon_e'].str.contains('NO_EXON_AFFECTED')]
        Structural_VCF = Structural_VCF[~Structural_VCF['Gene_name_g'].str.contains('NO_GENE_AFFECTED')]

        #Get the transcripts affected by these CNVs
        set_transcripts = set(Structural_VCF['Transcript_stable_ID_g'])

        #Find the location of all exons of these transcripts
        exons_data = import_data.loeuf_exons()
        exons_data = exons_data.loc[exons_data['Transcript_stable_ID'].isin(set_transcripts)]

        # 2 or more transcripts of the same gene can have common exons
        # we need to delete repeat exons to not call the same SNPs many times
        exons_data = exons_data.drop_duplicates(subset= ['Exon_region_start', 'Exon_region_end','Exon_rank_in_transcript'],
                                                keep = 'last').reset_index(drop = True)

        # Now we have the coordinates of all exons of affected genes
        # We can use this to filter SNPs that overlap with these coordinates
        # To do this, we need to make a bed file to do the filter
        output_path = '/re_gecip/machine_learning/LOEUF_classifier_tool/Data_output/tmp/exons_affected_mapping_SNPs.bed'
        exons_data.to_csv(output_path,mode='w', sep='\t', index=False)

        # Filter and take short variants overlapping the exons of our transcripts affected by CNVs
        # Subprocess needs to be dinamic because the -a parameter changes per run
        first_part = 'bedtools intersect -wb -header -b /re_gecip/machine_learning/LOEUF_classifier_tool/Data_output/tmp/exons_affected_mapping_SNPs.bed -a '
        third_part = '> /re_gecip/machine_learning/LOEUF_classifier_tool/Data_output/tmp/filtered_SNPs_1.vcf'
        to_subprocess = first_part + str(Standard_VCF_path) + third_part

        try:
            subprocess.run(to_subprocess, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(e.output)
            log.log_messages('ERROR:        Bedtools didnt work line 99 short_variant.py)')
            log.log_messages('INFO:         Annotating genes by alternative route')
            data_presentation2.hand_snp_annotation_error(genes,file_name)
            quit()



class Snp_annotation:
        def __init__(self,genes,file_name):
                self.genes = genes
                self.file_name = file_name
        def Annotate_snp(self, genes=None):
                '''
                Run a bash script to get high confident variants
                Annotate these variants to get transcript ID
                Merge this with structural data
                to add SNPs with affected transcripts
                '''
                genes=self.genes
                file_name = self.file_name
                log.log_messages('INFO:        Starting LOFTEEE_v99.1 annotation (this may take a while)')
                # VEP must be run in a bash script
                # VEP produce a lot of messages in the terminal. I have saved this in a file
                # with &>
                subprocess.run('/re_gecip/machine_learning/LOEUF_classifier_tool/Application/classes/SNP_annotation.sh &> /re_gecip/machine_learning/LOEUF_classifier_tool/Data_output/SNP_annotation.err', shell=True)
                # Load variants annotated into a data frame

                # Descard metadata
                # And ensure we always get the same columns
                columns2type ={'Uploaded_variation': str,
	                           'Location': str,
			                   'Allele': str,
			                   'Gene': str,
			                   'Feature': str,
			                   'Feature_type': str,
			                   'Consequence': str,
			                   'cDNA_position': str ,
			                   'CDS_position': str,
			                   'Protein_position': str ,
			                   'Amino_acids': str,
			                   'Codons': str,
			                   'Existing_variation': str,
			                   'Extra': str}

                log.log_messages('SUCCESS:     VEP annotation done')

                # HC_VEP_SNPs_annotated.vcf contains only High confident LoF variants
                # The number of HC LoF variants is low so some genomes could have 0 HC loF variants in ROI
                # We need to deal with this otherwise errors will occuers downstream

                hc_lof_snp = '/re_gecip/machine_learning/LOEUF_classifier_tool/Data_output/tmp/filtered_and_annotated_VEP_HC_SNPs.vcf'

                if os.stat(hc_lof_snp).st_size == 0:
                    # If size is = 0 that means there is not HC LoF variants in ROI
                    # So, inform the user about this and takes alternative route
                    log.log_messages('INFO:         Number of High Confident compound heterozigous variants: 0')
                    log.log_messages('INFO:         Annotating genes by alternative route')
                    data_presentation2.hand_snp_annotation_error(genes,file_name)

                    # stop the function here because the alternative route not only manages when there is not HC LoF variants
                    # but also does what the main route does such as the mane staff
                    return

                    # If size is > 0 at least one HC LoF variants has been found
                    # so, keep going with improvement data presentation 
                else:
                    snp_loftee = pd.read_csv(hc_lof_snp,
                                         sep = '\t',
                                         index_col=False,
                                         names = columns2type.keys(),
					 dtype = columns2type)

                    # Same variants may have more than one annotation.
                    # VEP output duplicate variants with more than one annotation
                    # Delete duplicate but leave the one that have different consequence
                    snp_loftee = snp_loftee.drop_duplicates(
                        subset = ['Location','Consequence'],
                        keep = 'last').reset_index(drop=True)
                    snp_loftee[['#Chro','Start']] = snp_loftee['Location'].str.split(':',expand = True)
                    snp_loftee['End'] = snp_loftee['Start']

                    # Get genotype
                    snp_loftee['GT'] = snp_loftee['Extra'].apply(lambda x: x.split(';')[1])
                    snp_loftee['HGVSc'] = snp_loftee['Extra'].apply(lambda x: x.split(';')[4])
                    
                    # Save in statistics number of High Confident compound heterozigous variants
                    hclof_variants = str(len(snp_loftee.HGVSc.unique()))         
                    log.fill_statistics_file(hclof_variants)
                    text = "INFO:        Number of High Confident compound heterozigous variants: {}".format(hclof_variants)
                    log.log_messages(str(text))

                    # Now, SNPs and genes affected are linked in the same data frame
                    snp_genes = genes.merge(snp_loftee,
                                        how='left',
                                        left_on='Gene_stable_ID_g',
                                        right_on='Gene')
                    


                # Requested by Clinical scientist to add MANE and gene-phenotype association
                # Instead of adding and changing the transcripts dataset (which will need revalidation)
                # I have created a independent dataset with the same transcripts ID attached Mane + phenotype
                # This is also a convinient method because if we increase in columns the size of the exons dataset
                # pandas would become slower
                # The mane df also contain Transcript length, Interpro start and Interpro end
                # These will be neccesary for sorting genes (will be explained later)

                # Load this dataset
                mane_file = '/re_gecip/machine_learning/LOEUF_classifier_tool/Resources/For_annotation/trancript_MANE_phenotype_v1.1.tsv'
                mane = pd.read_csv(mane_file,
                                         sep = '\t',
                                         index_col=False,
                                         dtype=str)
                # Merge by trancript ID
                snp_genes_mane = snp_genes.merge(mane,
                                             how = 'left',
                                             left_on = 'Transcript_stable_ID_g',
                                             right_on = 'Transcript stable ID')
                return snp_genes_mane
