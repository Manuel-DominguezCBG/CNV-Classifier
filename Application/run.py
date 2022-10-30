
'''
Name:            main.py
Author:          manue.dominguezbecerra@nhs.net
Description:     The highest script in this program.
                 The purpose of this script is to trigger the run.py script
                 and allow us to handle errors and exceptions. 
                 Also to carry out cleaning tasks whatever happens during the workflow
States:          RELEASE

Change control

CNV-Classifier v1.1.0-a.1
The program is able to do CNV annotation (genes and transcripts).

CNV-Classifier v1.1.0-a.2
The program is able to compound heterozygous annotation. Log file and statistics dataset features added.

CNV-Classifier v1.1.0-a.3
Main and complementary reports are now created as the user requested. Automating.sh allow the user to analyse >1 genome.

CNV-Classifier v1.1.0-rc.3
Minor changes before Release. Added The control-run test.

CNV-Classifier v1.1.0
First release. 

'''
__version__ = '1.1.0'

# Import libraries modules
import sys
import subprocess
from classes.main import CnvClassifier
from classes import log

try:

    # Create log file 
    log.log_messages('INFO:        CNV-Classifier_v1.1.0 STARTING')
    # Take input user
    user_input = sys.argv[1]
    # Run the CNV-Classifier
    run = CnvClassifier(user_input)
    run.run_application()
    log.log_messages('INFO:         Application finished correctly')
except Exception as e:
    # Catching too general exception but this is ok here
    # as many type of error can happens and we cannot hand all type of exception
    # at this tope level
    log.log_messages('ERROR:        Somenthing went wrong ')
    TEXT = f"ERROR:        {e}"  # Using f-string
    log.log_messages(str(TEXT))
finally:
    # do cleanup whatever happens
    print('Doing some cleanup')
    log.end_statisticts_file()
    subprocess.run('rm /re_gecip/machine_learning/LOEUF_classifier_tool/Data_output/tmp/*',
                   shell=True,
                   check=True
                   )


# This 3 lines also run the code and if a error occurs,
# this gives more details than the try/except block. 
# never used but in case future developers needes
# Only for development as no tmp files will be deleted (!)
'''
user_input = sys.argv[1]
run = CnvClassifier(user_input)
run.run_application()
'''
