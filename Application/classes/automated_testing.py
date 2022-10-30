'''
Name:           automated_testing.py
Author:         manuel.dominguezbecerra@nhs.net
Description:    The first approach we developed in this project
                to carry out annotation is this script. I tested that this
                really does what I want. However, I found that interset (bedtools)
                is a more efficient way to approach this. Therefore, I used bedtools
                and then I left this to validated the first bedtools way. 
                For 1-genome analysis, this delaes the execution a 0.68 extra secons what is ok
                When analysing a large number of genomes, this makes the full analysis very slow.
                However, during the development of this project, we found that we need to use HPC
                cluster which is a more powerful computer than the sourced asigned for the terminal
                of the remote linux system.
States:         RELEASE
'''

import os
import sys
import warnings
from classes import log


# Take final output and by an alternative procedure
# map CNVs and compare results
# Generate results in in a testing report
# Testing is incorporated in the application


class TestingAnnotation():
    '''
    To test in every execution 
    if annotation does what would expedted
    '''
    def __init__(self,cnvs):
        self.cnvs = cnvs

    def test_gene_annotation(self,cnvs=None):
        '''
        Is gene overlapping correct?
        This was actually a interset I developed before I found that
        interset (Bedtools) is faster and requiered
        less resources than this. So, I am reusing this code for validation
        '''
        cnvs=self.cnvs
        # Check if the overlapping between gene and CNV coordinates are correct, 
        # This return TRUE if correct. This is record in the CNVs_annotated_by_trancripts
        # so that we can see what is the problematic CNV (If any)
        cnvs['pass_validation_gene'] = (cnvs['START'].between(cnvs['Transcript_start_g'],
                                                               cnvs['Transcript_end_g']) | cnvs['END'].between(cnvs['Transcript_start_g'],
                                                               cnvs['Transcript_end_g'])) | (cnvs['Transcript_start_g'].between(cnvs['START'],
                                                               cnvs['END']) | cnvs['Transcript_end_g'].between(cnvs['START'],
                                                               cnvs['END']))
        # This approach return  FALSE when a CNV doesn't affect a gene
        # For these cases and only for these cases, convert this FALSE in TRUE
        cnvs.loc[cnvs.Transcript_stable_ID_g=='NO_GENE_AFFECTED', 'pass_validation_gene'] = True

        # If a FALSE is found, this means that the overlapping carry out in this test do not match
        # with the one carry out by Bedtools

        # Inform in log file the % of TRUEs found. 100% TRUE means, CNV mapping genes PASS validation.

        log.log_messages('VALIDATION:  Gene annotation ok? {}% True'.format(cnvs.pass_validation_gene.astype(str).value_counts(normalize=True).mul(100).round(1).astype(str).get("True")))

        return cnvs

    def test_exon_annotation(self,cnvs=None):
        '''
        Same approach than before but with exons
        '''
        cnvs=self.cnvs

        cnvs['pass_validation_exon'] =  (cnvs['START'].between(cnvs['Exon_region_start_e'],
                                                               cnvs['Exon_region_end_e']) | cnvs['END'].between(cnvs['Exon_region_start_e'],
                                                               cnvs['Exon_region_end_e'])) | (cnvs['Exon_region_start_e'].between(cnvs['START'],
                                                               cnvs['END']) | cnvs['Exon_region_end_e'].between(cnvs['START'],
                                                               cnvs['END']))

        cnvs.loc[cnvs.Constitutive_exon_e=='NO_EXON_AFFECTED', 'pass_validation_exon'] = True

        log.log_messages('VALIDATION:  Exon annotation ok? {}% True'.format(cnvs.pass_validation_exon.astype(str).value_counts(normalize=True).mul(100).round(1).astype(str).get("True")))

        return cnvs

def validation_data_loaded(data_loaded):
    '''
    Check columns are loaded in the correct order
    according to VCF specifications v.4.0, 4.1 and 4.2
    '''
    header_input = list(data_loaded.columns)
    header_input.pop(9) 
    if header_input != ['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']:
        log.log_messages('ERROR: Structure data_loaded failed')
        # If this is not correct very likely errors will appear. Better stop here 
        sys.exit(" Validation data loaded FAIL. App END")

def validation_structural_vcf(file_path):
    '''
    Structural vcf file when created, a "SV" is added 
    at the end of the file name (according to GEL)
    so check that we are taking a SV file
    This is very important when analysing several vcf files
    '''
    file_name = os.path.basename(file_path)
    if '.SV.vcf.gz' in 	file_name:
        log.log_messages('VALIDATION:  Input is a structural VCF file')
        # Add file name into log file
        message = 'INFO:        Sample running: ' + file_name
        log.log_messages(message)
    else:
        log.log_messages('VALIDATION-WARNING: Not a SV VCF file')

def filename_samplename(header_input, file_name):
    '''
    Validate sample name match with file name
    '''
    if file_name.split('.SV.',1)[0] == header_input[9]:
        log.log_messages('VALIDATION:  File name and sample name match')
    else:
        log.log_messages('VALIDATION-WARNING: File_name and sample_name doesnt match')
def columns_df(header_input):
    '''
    Validate that the column names of the created df are correct
    '''
    header_input.pop(9) 
    if header_input != ['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']:
        log.log_messages('ERROR: Structure data loaded failed')
        log.log_messages('ERROR: Application STOPS')
        sys.exit("It doesn't pass validation, check header VCF file")

def structural_df(vcf_df):
    '''
    Validated that the user input is a structural VCF file
    or at least that the file has been called appropriately
    '''
    if vcf_df.shape[1] != 10 and vcf_df.columns != []:
        log.log_messages(f" VALIDATION-WARNING: Unrecognised data frame shape in file: {file_name}")
        sys.exit("It doesn't pass validation 4 (read_file.py), check header VCF file")
        # VALIDATE1,2 don't stop the application because it doesn'f affect posterior analysis
        # However, validation 3 and 4 are critical and if header is not correct, it will produce error later
