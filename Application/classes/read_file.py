'''
Name:          read_file.py
Author:        manuel.dominguezbecerra@nhs.net
Description:   Convert VCF files into data frame
               and validate that data input and data output are correct.
States:        RELEASE
'''


# Import libraries
import io
import os
import gzip
import pandas as pd
from classes import log
from classes.automated_testing import *

class ReadFile():
    '''
    This class read a VCF file
    and does some data manipulation
    the outout is the full data found
    in the input of this class
    the filtering process happens
    in the following step
    '''
    def __init__(self,file_path):
        '''
        This is the built-in constructor method
        '''
        self.file_path = file_path

    def load_data(self):
        '''
        1) Convert VCF file into  data frame
           Read  header of the body dinamically and assign dtype
           This increase efficiency.
        2) Add some data to statistics
        3) Validate input data and the output of this function
        '''
        # Validate that input is a SV VCF file
        validation_structural_vcf(self.file_path)

        # Open the VCF file and read line by line
        with io.TextIOWrapper(gzip.open(self.file_path,'r')) as f:

            lines =[l for l in f if not l.startswith('##')]
            # Identify columns name line and save it into a dict
            # with values as dtype
            dinamic_header_as_key = []
            for liness in f:
                if liness.startswith("#CHROM"):
                    dinamic_header_as_key.append(liness)
                    # Declare dtypes
            values = [str,int,str,str,str,int,str,str,str,str]
            columns2detype = dict(zip(dinamic_header_as_key,values))

            vcf_df = pd.read_csv(
                io.StringIO(''.join(lines)),
                dtype=columns2detype,
                sep='\t'
            ).rename(columns={'#CHROM':'CHROM'})

            # Add sample name into statistics
            file_name = os.path.basename(self.file_path)
            log.fill_statistics_file(file_name)

            # Add number of variants into statistics
            log.fill_statistics_file(str(vcf_df.shape[0]))

            # Validate sample name match with file name:
            header_input = list(vcf_df.columns)
            filename_samplename(header_input, file_name)

            # Validate column names
            columns_df(header_input)

            # Validate vcf structure:
            structural_df(vcf_df)

            # Validate data loaded
            validation_data_loaded(vcf_df)
            return vcf_df
