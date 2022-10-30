'''
Name:            import_data.py
Author:          manuel.dominguezbecerra@nhs.net
Description:     A script to import requiered data
States:          INACTIVE Since bedtools is used, 
                 this script is not needed.
'''


import pandas as pd

def loeuf_genes():
    '''
     This function import the dataset we have created
	 with the LOEUF score and the genes from Ensembl.
    '''
    
    # Import LOEUF genes data set	    
    LOEUF_genes = pd.read_csv('/re_gecip/machine_learning/LOEUF_classifier_tool/Resources/For_annotation/LOEUF_GENES_COORDINATES_v2.tsv',
	                          sep = '\t',
							  names=['Chromosome',
							         'Transcript_start',
									 'Transcript_end',
									 'Transcript_stable_ID',
									 'canonical',
									 'Gene_stable_ID',
									 'Gene_name'])
    return LOEUF_genes


def loeuf_exons():
	'''
	
	'''
	    
	# Import LOEUF exons data set
	LOEUF_exons = pd.read_csv('/re_gecip/machine_learning/LOEUF_classifier_tool/Resources/For_annotation/LOEUF_EXONS_COORDINATES_v2.tsv',
	                          sep = '\t',
							  dtype={'#Chromosome':str,
							         'Exon_region_start':int,
									 'Exon_region_end':int,
									 'Exon_rank_in_transcript':int,
									 'Exon_stable_ID:':str,
									 'Constitutive_exon':int,
									 'Transcript_stable_ID':str,
									 'Gene_name':str,
									 'oe_lof_upper':str,
									 'oe_lof_upper_bin':str})
	return LOEUF_exons
'''
def loeuf_metrics_data():
	LOEUF_metrics_data = pd.read_csv('/re_gecip/machine_learning/LOEUF_classifier_tool/Resources/supplement/supplementary_dataset_11_full_constraint_metrics.tsv', sep = '\t', low_memory=False )
	
	# Take only what is needed so far
	# transcript to identify metrics
	# oe_lof_upper this is LOEUF score
	# Change this, dont limited the output, select what you need later
	return LOEUF_metrics_data[['transcript','oe_lof_upper']]
	# Data frame needs some changes
'''