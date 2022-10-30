'''
Name:            data_presentation2.py
Author:          manuel.dominguezbecerra@nhs.net
Description:     This is the alternative data_presentation script.
                 When the genome analysed does't present any high confident
                 predicted compound heterozigous, the pipeline returns an error
                 and the main report is not generated. 
                 To avoid this, this script achieves this task.
States:          RELEASE
'''

import os
import sys
import subprocess 
import pandas as pd
from classes.data_presentation import ImproventData
from classes import log

def hand_snp_annotation_error(genes,file_name):
    '''
    This function does exacly the same that
    the function found in data_presentation.py
    A separate function is needed because data in 
    both cases differ. For clarity, I have put
    this block of code in a separated script
    '''

    # Load the main and gene phenotype dataset
    mane_file = '/re_gecip/machine_learning/LOEUF_classifier_tool/Resources/For_annotation/trancript_MANE_phenotype_v1.1.tsv'
    mane = pd.read_csv(mane_file,
                        sep = '\t',
                        index_col=False,
                        dtype=str)
    # Merge by trancript ID
    genes = genes.merge(mane,
                            how = 'left',
                            left_on = 'Transcript_stable_ID_g',
                            right_on = 'Transcript stable ID')
    # List of columns to adapt 
    empty_columns = ['Uploaded_variation',
                    'Location',
                    'Allele',
                    'Gene',
                    'Feature',
                    'Feature_type',
                    'Consequence',
                    'cDNA_position',
                    'CDS_position',
                    'Protein_position',
                    'Amino_acids',
                    'Codons',
                    'Existing_variation',
                    'Extra',
                    '#Chro',
                    'Start',
                    'End',
                    'GT',
                    'HGVSc']
    genes[[empty_columns]] = ''
    # Add 0 to the number of HC LoF short variants
    log.fill_statistics_file('0 due to a error')
    output = ImproventData(genes, file_name)
    # Apply the data improvement as the main route
    output = output.reliable()
    # Add results to the file as the main route
    output_path = '/re_gecip/machine_learning/LOEUF_classifier_tool/Data_output/RESULTS/CNVs_annotated_by_genes.tsv'
    output.to_csv(output_path, sep='\t',mode='a', header=not os.path.exists(output_path))
