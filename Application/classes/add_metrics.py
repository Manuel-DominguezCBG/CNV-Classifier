'''
Name:            add_metrics.py
Author:          manuel.dominguezbecerra@nhs.net
Description:     Add LOEUF metrics
States:          RELEASE
'''
import pandas as pd
from classes import import_data


class LoeufMetrics():
    ''' 
    LOEUF metrics was already add when we map exons
    because the dataset I created has added LOEUF metrics.
    This dataset was created after I did this script because I found
    problems between the LOUEF dataset and the dataset used for annotation.
    Some gene names between versions were different.
    LOEUF score was added directly but I have reused this to validated that
    when I created the new dataset, all match between the new dataset used in annotation
    and the original LOEUF dataset
    '''
    def __init__(self,cnvs):
        self.cnvs = cnvs

    def oe_lof_upper(self,cnvs=None):
        cnvs=self.cnvs
        LOEUF_metrics_data = import_data.loeuf_metrics_data()
        cnvs = pd.merge(cnvs,LOEUF_metrics_data,left_on='Transcript_stable_ID_g',right_on='transcript', how='left')
        return cnvs
