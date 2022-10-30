'''
Name:            link.py
Author:          manue.dominguezbecerra@nhs.net
Description:     Tidy up data to link it with the 
		         following objetive.
                 From Sprint 3 to Sprint 4 the project was reorganised
                 this script link the workflow between sprint 3 (the cnv analysis)
                 with the Sprint 4 (the compound heterozigous analysis.)
                 This script also satisfy some user requeritments suggested during the development
States:          RELEASE
'''


import pandas as pd
import subprocess 
from classes import log
from classes import import_data
from classes import log



class OrganizateData():
    '''
    This scripts try to satisface 
    clinical users, deleting CNVs that 
    don't affect gen and/or exons.
    Then, annotate genes intead of CNVs
    '''
    def __init__(self,cnvs):
        self.cnvs = cnvs

    def cnv2genes(self,cnvs=None):
        '''
        '''
        cnvs=self.cnvs

        # Delete CNVs that affect ROI
        cnvs = cnvs[~cnvs['Constitutive_exon_e'].str.contains('NO_EXON_AFFECTED')]
        cnvs = cnvs[~cnvs['Gene_name_g'].str.contains('NO_GENE_AFFECTED')]
        cnvs['canonical_g'] = cnvs['canonical_g'].replace(['True','False'],['(Canonical)',''])
        # Most of the info generated is irrelevant for clinicians, let's take 
        # what the told me they want
        # I put in a list what they need
        clinicians_need = ['Gene_name_g',
                           'Gene_stable_ID_g',
                           'Transcript_stable_ID_g',
                           'canonical_g','CNV_TYPE',
                           'Exon_rank_in_transcript_e',
                           'oe_lof_upper_e',
                           'oe_lof_upper_bin_e',
                           'amount_overlap_g']

        cnvs = cnvs[clinicians_need]
        genes = cnvs.groupby(['Gene_name_g',
                              'Gene_stable_ID_g',
                              'Transcript_stable_ID_g',
                              'canonical_g',
                              'CNV_TYPE',
                              'oe_lof_upper_e',
                              'oe_lof_upper_bin_e',
                              'amount_overlap_g'
                             ])['Exon_rank_in_transcript_e'].apply(list).reset_index(name='Exons ranks affected')
        return genes



