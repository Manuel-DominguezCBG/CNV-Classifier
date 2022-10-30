'''
Name:            main.py
Author:          manuel.dominguezbecerra@nhs.net
Description:     This scripts takes the user input from run.py
                 If the user input is run_control, a validation process is executed
                 by running the function control_vcf()
                 If the user input us a SV VCF, then the main analysis is carried out
                 by running, ReadFile, Filter and so on.

States:          RELEASE
'''

# Import libraries modules
import os
from classes.read_file import ReadFile
from classes import log
from classes.filter import Filter
from classes.annotation import Annotation
from classes.add_metrics import LoeufMetrics
from classes.automated_testing import *
from classes.control_sample import *
from classes.link import OrganizateData
from classes.short_variants import *
from classes.data_presentation import ImproventData


class CnvClassifier:
    '''
    Call all functions and classes needed to run
    the application.
    '''

    def __init__(self,user_input):
        '''
        This is the built-in constructor method
        Save the user input (file path + SV VCF file name)
        '''
        self.user_input = user_input
        
    def run_application(self, user_input=None):
        '''
        This function run the application
        and the control sample
        '''
        user_input=self.user_input

        ########################################################
        #    Run test control:
        #                To run a manually created VCF file
        #                To ensure application does what expected
        #                Details in README file
        #                To run this >python run.py run_control<
        #
        ########################################################

        if user_input == 'run_control':

            # Run test control
            control_vcf()
            sys.exit()
            #quit()


        ########################################################
        # Standard pipeline   
        #    First part:
        #                Identify CNVs  from a VCF file
        #                Annotate genes and exons
        #                Add LOEUF metrics
        #                Organizate data (DEV)
        #
        ########################################################


        # Create statistics.csv file and variable builder
        log.write_statistics_file()

        # Generated a data frame from the VCF file
        data = ReadFile(user_input)
        data_loaded = data.load_data()

        # Filter and take CNVs
        cnv_df = Filter(data_loaded)
        cnv_df = cnv_df.filter_vcf_return_cnv()

        # Annotate genes and trancripts
        cnv_mapped_genes = Annotation(cnv_df)
        cnv_mapped_genes = cnv_mapped_genes.gene_annotation_bedtools()

        # Annotate exons
        cnv_mapped_exons = Annotation(cnv_mapped_genes)
        cnv_mapped_exons = cnv_mapped_exons.exon_annotation_bedtools()

        # Test  annotation genes
        # Same annotation done with a second method to test
        testing_genes_annotation = TestingAnnotation(cnv_mapped_genes)
        testing_genes_annotation = testing_genes_annotation.test_gene_annotation()

        # Test  annotation exons
        testing_exons_annotation = TestingAnnotation(cnv_mapped_exons)
        testing_exons_annotation = testing_exons_annotation.test_exon_annotation()

        # Add LOEUF metrics
        cnv_mapped_metrics = LoeufMetrics(cnv_mapped_exons)
        cnv_mapped_metrics = cnv_mapped_metrics.oe_lof_upper()

        # CNVs are annotated by genes in a tsv file
        # However, it was requested by the clinical scientist team
        # to provide full details about CNVs
        # such as CNVs coordinates of all affected exons.
        # This information is saved in CNVs_annotated_by_trancripts

        output_path = '/re_gecip/machine_learning/LOEUF_classifier_tool/Data_output/RESULTS/CNVs_annotated_by_trancripts.csv'
        cnv_mapped_metrics.to_csv(output_path,
                                  mode='a',
                                  header=not os.path.exists(output_path))


        # Organised CNVs by Gene
        genes = OrganizateData(cnv_mapped_metrics)
        genes = genes.cnv2genes()

        ########################################################
        #    Second part:
        #             Take the counterpart VCF file
        #             Identify SNPs and indels that overlap with transcripts of affected genes
        #             The compuond heterozigous analysis
        ########################################################

        log.log_messages('INFO:        Looking for Compound Heterozygous Variants')

        # Find in GEL database the counterpart VCF file (called standard VCF)
        # See README for explanation of how VCF files are organised per participant in GEL
        standard_vcf_path = user_input.replace(".SV","")
        # This is the absolute path of Standard VCF file

        # Identify SNPs and indels that overlap with affected transcripts
        # See README file for what was specifically requested by the users

        file_name = os.path.basename(user_input).split('.SV',1)[0]
        snp_select = Snp_manipulation(standard_vcf_path,cnv_mapped_metrics, genes, file_name)
        snp_select = snp_select.Identify_affected_transcripts()

        # Annotate SNPs and indels and match genes data with short variants

        snps_and_genes = Snp_annotation(genes,file_name)
        snps_and_genes = snps_and_genes.Annotate_snp()


        ########################################################
        #    Part 3:
        #             Improve data prersentation
        #             As requested by clinicians and clinical scientists
        #
        ########################################################

        # Improve data presentation

        output = ImproventData(snps_and_genes, file_name)
        output = output.reliable()

        # This creates the primary output
        # Readme.md explain this file this as well as all outputs developed by this app  
        output_path = '/re_gecip/machine_learning/LOEUF_classifier_tool/Data_output/RESULTS/CNVs_annotated_by_genes.tsv'
        # CNVs_annotated_by_genes is the main report of this project
        output.to_csv(output_path, sep='\t',mode='a', header=not os.path.exists(output_path))
