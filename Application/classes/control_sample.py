'''
Name:           control_sample.py
Author:         manuel.dominguezbecerra@nhs.net
Description:    This script is to carry out some validation procedure
                by the user in a easy way.
States:         RELEASE
'''

from classes.read_file import ReadFile
from classes import log 
from classes.filter import Filter
from classes.annotation import Annotation
from classes.add_metrics import LoeufMetrics
from classes.automated_testing import *
from classes.link import OrganizateData
from classes.short_variants import *
from classes.data_presentation import ImproventData

def control_vcf():
    ''' 
    Test by using a control VCF file
    2 manually VCF files have been created to check that the application 
    does what expected
    /re_gecip/machine_learning/LOEUF_classifier_tool/run_control

    The SV VCF file contains 'fake' CNVs. These have been created in specific coordinates
    to check that CNV annotation works in different scenarios

    Ths standard VCF file contains several 'fake' SNPs and indels located in specific coordinates
    However, two are called by LOFTEE and 1 is added to the final report
    This has been done deliberately
    to check that these variants are filtered correctly based on users' requeritments.

    (These are not provided in Github because the Standard VCF file also provides real variation, 
    and GEL would refuse them.)

    After creating these VCF files, we have run this application and
    the output has been checked manually.

    When the user type python run.py run_control,
    the application recognises this and runs this function which basically runs
    the full application and at the end checks that the output is correct.

    When I checked the output manually I created the md5 bashed 
    of the two main output files the main and the comp. reports.
    When the user runs this control test, this function generates new MD5 hashes
    and check that the one I created and the one generated during the execution are the same.

    This is an easy test the user can do to verify that everything is ok.

    This is important because the HPC computer used in GEL is in constant change
    and errors outside my control can happen.
    For example, during development, I have seen how the name of the modules needed here have changed
    (E.g. bio/BEDTools/2.19 changed to bio/BEDTools/2.19.1).
    If changes like this happen again and affect this program,
    the user can verify that at least the app. does what is expected.

    This method has the drawback that each new development that affects
    the final result requires updating the sample control. This increases the workload
    '''

    log.log_messages('--------VALIDATION STARTING: Running control VCF file---------')

    # Take user input (PATH/filename)
    # This is replaced by the control file
    user_input ='/re_gecip/machine_learning/LOEUF_classifier_tool/run_control/MD5_checking/Control-DNA_A09.SV.vcf.gz'

    # Followings lines are copy and pasted from main.py	

    ## Create statistics.csv file and variable builder
    log.write_statistics_file()

    ## Generated a data frame from the VCF file
    data = ReadFile(user_input)
    data_loaded = data.load_data()

    ## Filter and take CNVs 
    cnv_df = Filter(data_loaded)
    cnv_df = cnv_df.filter_vcf_return_cnv()


    ## Annotate genes 
    cnv_mapped_genes = Annotation(cnv_df)
    cnv_mapped_genes = cnv_mapped_genes.gene_annotation_bedtools()

    ## Annotate exons 
    cnv_mapped_exons = Annotation(cnv_mapped_genes)
    cnv_mapped_exons = cnv_mapped_exons.exon_annotation_bedtools()

    ## Test  annotation genes
    testing_genes_annotation = TestingAnnotation(cnv_mapped_genes)
    testing_genes_annotation = testing_genes_annotation.test_gene_annotation()

    ## Test  annotation exons
    testing_exons_annotation = TestingAnnotation(cnv_mapped_exons)
    testing_exons_annotation = testing_exons_annotation.test_exon_annotation()

    ## Add LOEUF metrics
    cnv_mapped_metrics = LoeufMetrics(cnv_mapped_exons)
    cnv_mapped_metrics = cnv_mapped_metrics.oe_lof_upper()

    # Generate the file with the CNV annotated
    output_path = '/re_gecip/machine_learning/LOEUF_classifier_tool/Data_output/RESULTS/CNVs_annotated_by_trancripts.csv'
    cnv_mapped_metrics.to_csv(output_path,mode='a', header=not os.path.exists(output_path))

    # Now we continue with SNPs and indels annotation

    # Organised CNVs by Gene
    genes = OrganizateData(cnv_mapped_metrics)
    genes = genes.cnv2genes()

    Standard_VCF_path = user_input.replace(".SV","")
    # This is the absolute path of Standard VCF file

    file_name = os.path.basename(user_input).split('.SV',1)[0]
    snp_select = Snp_manipulation(Standard_VCF_path,cnv_mapped_metrics, genes, file_name)
    snp_select = snp_select.Identify_affected_transcripts()

    snps_and_genes = Snp_annotation(genes,file_name)
    snps_and_genes = snps_and_genes.Annotate_snp()

    # Improve data presentation

    output = ImproventData(snps_and_genes, file_name)
    output = output.reliable()

    # This creates the primary output
    # Readme.md file this as well as all outputs developed by this app  
    output_path = '/re_gecip/machine_learning/LOEUF_classifier_tool/Data_output/RESULTS/CNVs_annotated_by_genes.tsv'
    output.to_csv(output_path, sep='\t',mode='a', header=not os.path.exists(output_path))


    # Here the application end. 

    # Now, we generate the md5 bash of the 2 main outputs
    # and compare with the MD5 hashes I have created after checking manually that everything looks ok

    # The MD5 hashes created during development were
    known_hash_CNVs_annotated_by_genes = '7184bdc833eae8d550a25b6884593810'
    known_hash_CNVs_annotated_by_trancripts = '5f7e9de3eb9b932e2ec8d73d6cc4e6a8'

    # and I get the hash from the files generated during this control process
    
    generated_hash_CNVs_annotated_by_genes = subprocess.check_output('md5sum /re_gecip/machine_learning/LOEUF_classifier_tool/Data_output/RESULTS/CNVs_annotated_by_genes.tsv', shell=True).decode("utf-8").split(' ',1)[0]

    generated_hash_CNVs_annotated_by_trancripts = subprocess.check_output('md5sum /re_gecip/machine_learning/LOEUF_classifier_tool/Data_output/RESULTS/CNVs_annotated_by_trancripts.csv', shell=True).decode("utf-8").split(' ',1)[0]

    # Finally, we heck them.

    checking = [known_hash_CNVs_annotated_by_trancripts == generated_hash_CNVs_annotated_by_trancripts,
                known_hash_CNVs_annotated_by_genes == generated_hash_CNVs_annotated_by_genes
                ]

    if all(checking):
        log.log_messages('--------VALIDATION-SUCCESS: md5 checking OK ---------')
    else:
        log.log_messages('--------VALIDATION-WARNING: md5 checking FAIL ---------')
    
    log.log_messages('-------- VALIDATION FINISHED ---------')
    
    # Remove the files created to not interfere with following runs
    subprocess.run('rm  /re_gecip/machine_learning/LOEUF_classifier_tool/Data_output/log.txt /re_gecip/machine_learning/LOEUF_classifier_tool/Data_output/statistics.csv', shell=True)
    subprocess.run('rm  /re_gecip/machine_learning/LOEUF_classifier_tool/Data_output/RESULTS/CNVs_annotated_by_trancripts.csv /re_gecip/machine_learning/LOEUF_classifier_tool/Data_output/RESULTS/CNVs_annotated_by_genes.tsv', shell=True)


'''
This is the output you should get if this run control test goes ok

2022-09-28 19:46:00,992 INFO:        CNV-Classifier_v1.3 STARTING
2022-09-28 19:46:00,992 --------VALIDATION STARTING: Running control VCF file---------
2022-09-28 19:46:01,017 INFO:        statistics.csv created
2022-09-28 19:46:01,017 VALIDATION:  Input is a structural VCF file
2022-09-28 19:46:01,017 INFO:        Sample running: Control-DNA_A09.SV.vcf.gz
2022-09-28 19:46:01,033 VALIDATION:  File name and sample name match
2022-09-28 19:46:01,047 SUCCESS:     CNVs filtered
2022-09-28 19:46:01,047 INFO:        Number of CNVs found in this VCF file are: 5
2022-09-28 19:46:01,474 SUCCESS:     CNV mapped (genes)
2022-09-28 19:46:03,460 SUCCESS:     CNV mapped (exons)
2022-09-28 19:46:03,464 VALIDATION:  Gene annotation ok? 100.0% True
2022-09-28 19:46:03,467 VALIDATION:  Exon annotation ok? 100.0% True
2022-09-28 19:46:04,814 INFO:        Standard VCF file exists
2022-09-28 19:46:05,901 INFO:        Starting LOFTEEE_v99.1 annotation (this may take a while)
2022-09-28 19:46:32,861 SUCCESS:     VEP annotation done
2022-09-28 19:46:32,885 INFO:        Number of High Confident compound heterozigous variants: 2
2022-09-28 19:46:33,150 --------VALIDATION-SUCCESS: md5 checking OK ---------
2022-09-28 19:46:33,150 -------- VALIDATION FINISHED ---------
Doing some cleanup
'''