'''
Name:            data_presentation.py
Author:          manuel.dominguezbecerra@nhs.net
Description:     Tidy up data to link it with the 
		         Improve how data is presented
                 following clinicians suggestions 
States:          RELEASE
'''
import pandas as pd
from classes import log
class ImproventData():
    '''

    '''
    def __init__(self,snps_and_genes,file_name):
        self.snps_and_genes = snps_and_genes
        self.file_name = file_name

    def reliable(self,snps_and_genes=None,file_name=None):
        '''
        Show data as requiered by clinicians
        and clinical scientists
        '''
        snps_and_genes = self.snps_and_genes
        file_name = self.file_name

        # There are repeated rows due to more than one phenotype per gene
        # Concatenate phenotype per gene and drop repeat columns

        # Canonical attached to transcript
        snps_and_genes['Transcript ID'] = snps_and_genes['Transcript_stable_ID_g'] + snps_and_genes['canonical_g']
        #snps_and_genes = snps_and_genes.drop(['Transcript_stable_ID_g','canonical_g'],axis=1)

        # Attach gene name to gene ID
        snps_and_genes['Gene'] = snps_and_genes['Gene_name_g'] + '--' + snps_and_genes['Gene_stable_ID_g']
        snps_and_genes = snps_and_genes.drop(['Gene_name_g','Gene_stable_ID_g'],axis=1)

        # Sort exon rank
        snps_and_genes['Exons ranks affected'] = snps_and_genes['Exons ranks affected'].apply(lambda lst: sorted(lst,key=int))
        snps_and_genes['Exons ranks affected'] = [','.join(map(str,l)) for l in snps_and_genes['Exons ranks affected']]

        # Sort by LOEUF metrics. More intolerant with SNPs/INDELS first
        # Take genes with short variats
        genes_with_snp = snps_and_genes[~snps_and_genes['HGVSc'].isnull()]
        snps_and_genes = snps_and_genes[snps_and_genes['HGVSc'].isnull()]
        
        # Sort them
        genes_with_snp = genes_with_snp.sort_values(by = 'oe_lof_upper_e',
                                                    ascending = False)
        # If any, take genes without LOUEF metrics
        genes_with_na = snps_and_genes[snps_and_genes['oe_lof_upper_e'].str.contains('Nan')]
        snps_and_genes = snps_and_genes[~snps_and_genes['oe_lof_upper_e'].str.contains('Nan')]

        # Sort the rest
        snps_and_genes = snps_and_genes.sort_values(by = 'oe_lof_upper_e',
                                                    ascending = False)

        # Concat, genes with short variants firts, then genes without SNPs and finally na genes
        snps_and_genes = pd.concat([genes_with_snp,snps_and_genes,genes_with_na])

        # Save 'Genes_affected' in the statistic file
        log.fill_statistics_file(str(len(list(snps_and_genes.Gene.unique()))))

        # Clinical scientis wants to know if mane and the canonical transcript
        # of our genes are affected. But also, they prefer one row per gene, 
        # Solution:
        '''
        From this

        Gene    Canonical     MANE
        Gene1                NA_0001
        Gene1   ENST00001    
        Gene2   
        Gene3                NA_0004

        To this

        Gene    Canonical     MANE
        Gene1     Yes          Yes
        Gene2      No          No
        Gene3      No          Yes
        '''
        # It is inevitable to reduce the amount of info displayed
        # when we improve presentation. But Clinical scientists and clinicians
        # are not interested on transcript info, so that is ok
  
        important_transcrips = lambda x: "Yes" if x.any() else "No" 

        # Genes are repeated if a gen has >1 phenotypes
        # To avoid repetition, we can concatenate phenotypes in only one row 
        phenotype = lambda x: ", ".join(y for y in x.fillna("").unique()if y)

        # We have been treating some numeric values as str to avoid error during this workflow
        # Now, I need to do some maths with them, we convert these to numeric

        numeric_columns = ['oe_lof_upper_e','amount_overlap_g']
        snps_and_genes[numeric_columns] = snps_and_genes[numeric_columns].apply(pd.to_numeric, errors = 'coerce', axis = 1)

        # Now, we can apply all and create a 1-row-per-gene data frame
        snps_and_genes = (snps_and_genes.groupby(["Gene" ]).agg({"CNV_TYPE": 'first',
                                                                 'canonical_g': important_transcrips,
                                                                 'RefSeq match transcript (MANE Select)': important_transcrips,
                                                                 'oe_lof_upper_e':['max','min'],
                                                                 'HGVSc':'first',
                                                                 'GT': 'first',
                                                                 'Consequence': 'first',
                                                                 'Uploaded_variation':'first' ,
                                                                 'Allele': 'first',
                                                                 'Phenotype description': phenotype,
                                                                 'amount_overlap_g':'max',
                                                                }))

        # Rename columns
        snps_and_genes.columns = ['CNV_TYPE',
                                  'canonical_g',
                                  'RefSeq match transcript (MANE Select)',
                                  'LOEUF metrics (max)',
                                  'LOEUF metrics (min)','HGVSc',
                                  'GT','Consequence',
                                  'Uploaded_variation',
                                  'Allele',
                                  'Phenotype description',
                                  'CNV overlaping with the gene']
        
        # Now, for presentation, revert from multiindex to single index data frame

        snps_and_genes = snps_and_genes.reset_index()

        # Clinical scientists suggested that if >1 genomes are analysed
        # intead of creating one file per genome, save all in one file
        # but adding identification column as index.
        snps_and_genes['Participant identification'] = file_name
        snps_and_genes = snps_and_genes.set_index('Participant identification')

        # An finally, sort affected genes by LOEUF metric (intolerant first)
        snps_and_genes = snps_and_genes.sort_values(by=['LOEUF metrics (max)'],ascending=True)


        return snps_and_genes

