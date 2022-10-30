'''
Name:            annotation.py
Author:          manue.dominguezbecerra@nhs.net
Description:     CNV annotation. This add gene and exons annotation
                 but the dataset we have created for this also provided
                 additional data (e.g. exon are linked with their transcripts
                 so this also provides transacripts annotation)
States:          RELEASE
'''

import subprocess
import warnings
import pandas as pd
from classes import log

class Annotation():
    '''
    Carry out cnv annotation via
	bedtools intersect.
	First function annotate genes
	second annotate exons and transcripts
	This class also record data to store in the statistics file
    '''
    def __init__(self,cnvs):
        self.cnvs = cnvs

    def gene_annotation_bedtools(self,cnvs=None):
        '''
        Map CNVs with genes
        This is done with interset (bedtools)
        '''
        cnvs=self.cnvs
        # Generate the bed file
        cnvs[['#CHROM','START','END','CNV_TYPE', 'CNV_ID']].to_csv('/re_gecip/machine_learning/LOEUF_classifier_tool/Data_output/tmp/cnvs_to_annotate_genes.bed',
              sep='\t',
              index=False)
        # Do bedtools intersect
        subprocess.run('bedtools intersect -wao -b /re_gecip/machine_learning/LOEUF_classifier_tool/Resources/For_annotation/LOEUF_GENES_COORDINATES_v2.tsv -a /re_gecip/machine_learning/LOEUF_classifier_tool/Data_output/tmp/cnvs_to_annotate_genes.bed > /re_gecip/machine_learning/LOEUF_classifier_tool/Data_output/tmp/cnv_mapped_genes.bed',
		               shell=True,
					   check=True)
        # Save results in a data frame
        cnvs = pd.read_csv('/re_gecip/machine_learning/LOEUF_classifier_tool/Data_output/tmp/cnv_mapped_genes.bed',
                           sep = '\t',
                           index_col=False,
                           names=['#CHROM',
                                  'START',
                                  'END',
                                  'CNV_TYPE',
                                  'CNV_ID',
                                  'Chromosome_g',
                                  'Transcript_start_g',
                                  'Transcript_end_g',
                                  'Transcript_stable_ID_g',
                                  'canonical_g',
                                  'Gene_stable_ID_g',
                                  'Gene_name_g',
                                  'amount_overlap_g'])

        # Interset adds "." and '-1' when CNVs do not overlap with any gene
        # Change this for something more informatived
        # Change this for No gene affected
        cnvs['Gene_name_g'] = cnvs['Gene_name_g'].str.replace('.','NO_GENE_AFFECTED')
        cnvs['Chromosome_g'] = cnvs['Chromosome_g'].str.replace('.','NO_GENE_AFFECTED')
        cnvs['Gene_stable_ID_g'] = cnvs['Gene_stable_ID_g'].str.replace('.','NO_GENE_AFFECTED')
        cnvs['Transcript_stable_ID_g'] = cnvs['Transcript_stable_ID_g'].str.replace('.','NO_GENE_AFFECTED')

        # Save data in statistics file
        # Number of CNVs that didn't affected genes
        tmp = cnvs.loc[cnvs['Gene_name_g'] == 'NO_GENE_AFFECTED']
        log.fill_statistics_file(str(len(tmp)))

        # Number of non-affected-genes CNVs (GAIN and LOSS)
        log.fill_statistics_file(str(tmp["CNV_TYPE"].value_counts().get("GAIN")))
        log.fill_statistics_file(str(tmp["CNV_TYPE"].value_counts().get("LOSS")))
        log.log_messages('SUCCESS:     CNV mapped (genes)')

        # How many trancripts affected
        filter = cnvs[~cnvs['Chromosome_g'].str.contains('NO_GENE_AFFECTED')]
        log.fill_statistics_file(str(len(filter)))

        # Count how many CNVs affeted one gene
        one = filter.groupby('CNV_ID')['Gene_stable_ID_g'].count().reset_index()
        one = one[one['Gene_stable_ID_g'] == 1].shape[0]
        log.fill_statistics_file(str(one))

        # Count how many CNVs affeted 2 genes
        two = filter.groupby('CNV_ID')['Gene_stable_ID_g'].count().reset_index()
        two = two[two['Gene_stable_ID_g'] == 2].shape[0]
        log.fill_statistics_file(str(two))

        # Count how many CNVs affeted 3 genes
        three = filter.groupby('CNV_ID')['Gene_stable_ID_g'].count().reset_index()
        three = three[three['Gene_stable_ID_g'] == 3].shape[0]
        log.fill_statistics_file(str(three))

        # Count how many CNVs affeted 4 genes
        four = filter.groupby('CNV_ID')['Gene_stable_ID_g'].count().reset_index()
        four = four[four['Gene_stable_ID_g'] == 3].shape[0]
        log.fill_statistics_file(str(four))

        # Count how many CNVs affeted 5 genes
        five = filter.groupby('CNV_ID')['Gene_stable_ID_g'].count().reset_index()
        five = five[five['Gene_stable_ID_g'] == 5].shape[0]
        log.fill_statistics_file(str(five))

        # Count how many CNVs affeted >5 genes
        sixmore = filter.groupby('CNV_ID')['Gene_stable_ID_g'].count().reset_index()
        sixmore = sixmore[sixmore['Gene_stable_ID_g'] >5].shape[0]
        log.fill_statistics_file(str(sixmore))

        # How likely more than 1 transcripts is affected when a CNVs affect a gene
        trancripts = filter.groupby('Gene_stable_ID_g')['Transcript_stable_ID_g'].count().reset_index()
        transcripts = (trancripts[trancripts['Transcript_stable_ID_g'] >1].shape[0]/trancripts.shape[0])*100
        log.fill_statistics_file(str(transcripts))

        return cnvs


    def exon_annotation_bedtools(self,cnvs=None):
        '''
        Map CNVs with exons
        This is done with interset (bedtools)
        '''
        cnvs=self.cnvs
        cnvs.to_csv('/re_gecip/machine_learning/LOEUF_classifier_tool/Data_output/tmp/cnvs_to_annotate_exons.bed',
		            header=True,
					sep='\t',
					index=False)
        subprocess.run('bedtools intersect -wao -b /re_gecip/machine_learning/LOEUF_classifier_tool/Resources/For_annotation/LOEUF_EXONS_COORDINATES_v2.tsv -a /re_gecip/machine_learning/LOEUF_classifier_tool/Data_output/tmp/cnvs_to_annotate_exons.bed > /re_gecip/machine_learning/LOEUF_classifier_tool/Data_output/tmp/cnv_mapped_exons.bed',
                       shell=True)
        # Save results in a data frame
        # As we know the dtype of this data, we can specify the dtype in advances with this dict.
        columns2type ={'#CHROM': str,
                       'START': str,
                       'END': str,
                       'CNV_TYPE': str,
                       'CNV_ID': str,
                       'Chromosome_g': str,
                       'Transcript_start_g': str,
                       'Transcript_end_g': str ,
                       'Transcript_stable_ID_g': str,
                       'canonical_g': str,
                       'Gene_stable_ID_g': str ,
                       'Gene_name_g': str,
                       'amount_overlap_g': str,
                       'Chromosome_e': str,
                       'Exon_region_start_e': str,
                       'Exon_region_end_e': str,
                       'Exon_rank_in_transcript_e': str ,
                       'Exon_stable_ID_e': str,
                       'Constitutive_exon_e': str,
                       'Transcript_stable_ID_e': str,
                       'Gene_name_e': str,
                       'oe_lof_upper_e': str,
                       'oe_lof_upper_bin_e': str,
                       'amount_overlap_e':  str}
        cnvs = pd.read_csv('/re_gecip/machine_learning/LOEUF_classifier_tool/Data_output/tmp/cnv_mapped_exons.bed',
                        sep = '\t',
                        index_col=False,
                        names = columns2type.keys(),
                        dtype = columns2type,
                     	usecols = range(23))

        cnvs['Transcript_stable_ID_e'] = cnvs['Transcript_stable_ID_g'].where(cnvs['Transcript_stable_ID_e'] .eq('.')).fillna(cnvs['Transcript_stable_ID_e'])
        '''
        Interset map canonical and no caninical exons twice
        E.g. GENE-1 has three exons, canonical transcripts are exons 1,2 and 3. 
		Exon 1 and 3 are also non-canonical exons.
        If  a CNV intersepts these 3 exons, the result would be.
        CNVID	Gene_name	Transcript	Canonical 	    Exon	Transcript
        CNV1	GENE-1		ENST0001	True			1		ENST0001
        CNV1	GENE-1		ENST0001	True			2		ENST0001
        CNV1	GENE-1		ENST0001	True			3		ENST0001
        CNV1	GENE-1		ENST0002	False		    1		ENST0002
        CNV1	GENE-1		ENST0002	False		    3		ENST0002
        CNV1	GENE-1		ENST0001	True			1		ENST0002
        CNV1	GENE-1		ENST0001	True			3		ENST0002
        The last 2 lines is info we dont need because we already know this
        To delete this "duplications", we just take the rows which Transcript_stable_ID_g == Transcript_stable_ID_e
        This will delete the ones which transcripts is not the same
        '''
        cnvs = cnvs[cnvs['Transcript_stable_ID_g'] == cnvs['Transcript_stable_ID_e']]

        # Some tidy up
        cnvs['Constitutive_exon_e'] = cnvs['Constitutive_exon_e'].str.replace('.','NO_EXON_AFFECTED')

        log.log_messages('SUCCESS:     CNV mapped (exons)')

        return cnvs
