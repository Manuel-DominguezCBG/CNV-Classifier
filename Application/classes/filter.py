'''
Name:          Filter.py
Author:        manue.dominguezbecerra@nhs.net
Description:   A script to select CNVs from the pool of structural variation.
               This also discards CNVs with low quality.
States:        RELEASE
'''
import sys
import configparser
import pandas as pd
from classes import log

# Read 'config.ini' file
config = configparser.ConfigParser()
config.read('/re_gecip/machine_learning/LOEUF_classifier_tool/Application/config.ini')
config.get('CNV_FILTER','QUAL')
config.get('CNV_FILTER','FILTER')
pd.options.mode.chained_assignment = None

class Filter:
    '''
    Take a data frame and select CNVs.
    Filter also by quality.
    Do some data manipulation
    Close application if not CNVs found (unlikely in GEL data)
    '''

    def __init__ (self,vcf_df):
        self.vcf_df = vcf_df

    def filter_vcf_return_cnv(self, vcf_df=None):
        '''
        Take only CNVs variation from the SV VCF file.
        This needs to be the first step because
        each type of SV presents different structure info
        and SV VCF are large files, therefore for a efficiency
        point of view this is also important
        '''
        vcf_df = self.vcf_df
        
        # This select CNVs only
        vcf_df['SVTYPE'] = vcf_df['INFO'].apply(lambda x: x.split(';')[0])
        cnv_df = vcf_df[vcf_df.SVTYPE == "SVTYPE=CNV"]
        # This action has a great inpact in this application and this approach
        # is explained in the documentation

        # Users requier to filter CNVs based on QUAL and FILTER
        # There is a configuration file in which user can change this values
        # Without the need to modify the code

        # This select variants which QUAL is > than QUAL
        cnv_df = cnv_df[cnv_df.QUAL > config.getint('CNV_FILTER','QUAL')]

        # This select variants that match with FILTER only!!
        # Variants that don't get the PASS flag have a very low QUAL
        # and we are not interested on variants with such a low QUAL.

        cnv_df = cnv_df[cnv_df.FILTER == config.get('CNV_FILTER','FILTER')]

        # Important values such END or CNV type are in the same column of the VCF
        # We split them  here
        cnv_df["CNV_TYPE"] = cnv_df["ID"].apply(lambda x:x.split(":")[1])
        cnv_df["START"] = cnv_df["ID"].apply(lambda x:x.split(":")[3])
        cnv_df["END"] = cnv_df["ID"].apply(lambda x:x.split(":")[4])

        # Add number of CNV into statistics
        log.fill_statistics_file(str(cnv_df.shape[0]))

        # Add number of gain and loss cnvs
        log.fill_statistics_file(str(cnv_df["CNV_TYPE"].value_counts().get("LOSS")))
        log.fill_statistics_file(str(cnv_df["CNV_TYPE"].value_counts().get("GAIN")))

        # As this will be the number of variants that will be analised.
        # Let's create a variant ID that take sample all important info we may need
	    # in posterior analysis for example

        pos = 9 # Column 9th is the name of the VCF ID
        colname = cnv_df.columns[pos]
        cnv_df = cnv_df.rename(columns={'CHROM':'#CHROM'})
        cnv_df['CNV_ID'] = colname +'|' + cnv_df['ID'].astype(str) +'|' + cnv_df['QUAL'].astype(str) +'|' + cnv_df['FILTER'] +'|' + cnv_df['FORMAT'] +'|' + cnv_df[cnv_df.columns[pos]]

	    # The hash is add at the beginning because after this,
        # a bed file will be created and comments in this file need
        log.log_messages('SUCCESS:     CNVs filtered')
        cnvs_number = len(cnv_df)
        if cnvs_number >=1:
            text = f"INFO:        Number of CNVs found in this VCF file are: {cnvs_number}"
            log.log_messages(str(text))
        else:
            # This is very unlikely, I haven't seen a VCF with less than 200 CNVs
            # But just in case to avoid this error
            log.log_messages('INFO: No CNVs found in this VCF file. Application finished')
            sys.exit('No CNVs found in this VCF file. Application finished')

        # User requiritments (2nd meeting): Put QUAL in the final output
        # Instead of adding a new column here which could need a lot of work to
        # accomodate that change across the code. Merge QUAL with CNV_TYPE
        # Then, at the end of the workflow, I split thet column to put QUAL in a new column

        cnv_df["CNV_TYPE"] = cnv_df["CNV_TYPE"] + "-" + (cnv_df['QUAL'].astype(str))
        return cnv_df[['CNV_ID','#CHROM','REF','ALT','QUAL','FILTER','CNV_TYPE','START','END',]]
