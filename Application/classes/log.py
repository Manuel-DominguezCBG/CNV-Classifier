'''
Name:          log.py
Author:        manuel.dominguezbecerra@nhs.net
Description:   This scripts achieve two aims.
               First, records ERROR, WARNING and INFO messages
               to explain if the app is working ok
               Second, creates the statistis file to record 
               genomics data    
States:        RELEASE
'''

import warnings
import logging
import csv
import subprocess
import os
import sys

log_file = '/re_gecip/machine_learning/LOEUF_classifier_tool/Data_output/log.txt'
statistics = '/re_gecip/machine_learning/LOEUF_classifier_tool/Data_output/statistics.csv'
tmp_confi = '/re_gecip/machine_learning/LOEUF_classifier_tool/Data_output/tmp/tmp_confi.txt' 
#### Using looging, a function can be done here to log any message in a log file

# To create and config logger
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(message)s',
    handlers=[
        logging.FileHandler(log_file),
        logging.StreamHandler(sys.stdout)])

# To create the object
logger = logging.getLogger()

# To set the threshold of logger to DEBUG
logger.setLevel(logging.DEBUG)

def log_messages(text):
    '''
    To be called in other scripts to record
    INFO/WARNING/ERROR messages 
    '''
    logging.info(text)

#### Statistics data

# It would be convinient to get data from the VCF files 
# when automatization the analysis of many VCF files
# as requested
def write_statistics_file():
    '''
    Create a csv file (if this doesn't exist)
    then add the header only when created
    This print messages that are in the code
    but if error occurs, this is reported in log.txt too.
    '''
    if os.path.exists(statistics) == False:
        try:
            with open(statistics, 'w') as f_object:
                writer_object = csv.writer(f_object)
                # The header of the statistics
                # If you would like to add more column, add the header here
                # the number of header and how many times this function is called
                # needs to the same, otherwise a error will occurs	
                header = ['File_name',
                          'Number_of_variants_in_SV-VCF',
                          'Number_of_CNV_in_SV-VCF',
                          'LOSS_CNV_per_vcf',
                          'GAIN_CNV_per_vcf',
                          'CNVS_not_affecting_genes',
                          'CNVS_not_affecting_genes_GAIN',
                          'CNVS_not_affecting_genes_LOSS',
                          'Number_of_transcripts_per_file',
                          'Number_CNVs_affecting_1_gene',
                          'Number_CNVs_affecting_2_genes',
                          'Number_CNVs_affecting_3_genes',
                          'Number_CNVs_affecting_4_genes',
                          'Number_CNVs_affecting_5_genes',
                          'Number_CNVs_affecting_more_than_5',
                          'How_likely_CNVs_affect_>1_transcripts_per_gene_affected_(%)',
                          'Number of High_confident_LoF_variants',
                          'Genes_affected',
                          ]
                writer_object.writerow(header)
                os.system('sed "$d" /re_gecip/machine_learning/LOEUF_classifier_tool/Data_output/statistics.csv')
        except OSError:
            logging.info("ERROR: Failed creating statistics.csv")
        else:
            logging.info("INFO:        statistics.csv created")


def fill_statistics_file(value):
    '''To record statistics data in the file'''
    with open(statistics, 'a') as f_object:
        f_object.write(value)
        f_object.write(',')

def end_statisticts_file():
    '''
    Add end of line when workflow is finised
    Data from the following VCF file (if any)
    will be recorded in a new line. 
    This way what we create is a data frame
    ready to be used. There are python libraires
    to create datasets but I haven't found in GEL
    '''
    with open(statistics, 'a') as f_object:
        f_object.write('\n')