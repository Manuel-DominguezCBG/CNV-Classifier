<<COMMENT
Name:            automating.sh
Author:          manuel.dominguezbecerra@nhs.net
Description:     This is an auxiliary script developed
                 to run the CNV-classifier over several genomes 
		         This script uses 'bsub' to submit a job for batched execution
                 by the lsbatch system.
                 Put the path/file_name of the genomes you wish to analysed 
                 in a file 
States:          RELEASE
COMMENT


##NECESSARY JOB SPECIFICATIONS
#BSUB -q long
#BSUB -L /bin/bash
#BSUB -P re_gecip_machine_learning
#BSUB -o '/re_gecip/machine_learning/LOEUF_classifier_tool/Data_output/job.%J.out'
#BSUB -e '/re_gecip/machine_learning/LOEUF_classifier_tool/Data_output/job.%J.err'
#BSUB -J JobName1
#BSUB -R "rusage[mem=100000] "
#BSUB -M 1
#BSUB -n 1
#BSUB -cwd your dir
#BSUB -B -N

# bsub documentation:
# https://www.ibm.com/docs/en/spectrum-lsf/10.1.0?topic=SSWRJV_10.1.0/lsf_command_ref/bsub.man_top.1.html

# To submit a job type 'bsub < automating.sh'  
# in the directory where automating.sh is

# To check if your application is running type 'bjobs -u your_user_name'

## LOAD MODULES AND CONDA ENV.
module load  bio/BEDTools/2.19.1
module load bio/BCFtools/1.9-foss-2019b

# GEL changes module names or locations
# (Happeed once during my development)
# If the error message comes up indicate that the module loaded doesn't exist
# Try 'module avail' and fing out what they did with the two modules needed
# or report an issue/raise a reques on 
# https://jiraservicedesk.extge.co.uk/plugings/servlet/desk

# Follow instrunction found in README.md file to generate the env. I have develop for this
# When done, activate it with these (change cnv by the name you asigned to your conda env.)
# README.md explains how to clone the cnv env. using requirements.txt
# This is to run ananconda
source /resources/conda/miniconda3/bin/activate
# This is to run the conda env created for this project
conda activate cnv

# This is a while loop to run CNV-Classifier as many times as genomes
# found in the file.
while read line  
do  
    python /re_gecip/machine_learning/LOEUF_classifier_tool/Application/run.py "$line" 
    # To know if all genomes found in the file have been analised
    COUNTER=$((COUNTER + 1))
    
done <  '/YOUR_PATH/file_name'

# Give permission to files created
# Both files are created in the first iteration
# Sometimes, when python tries to add data in the second iteration
# it has not permission. 700 protect the files agaist any access from other users
# while the issuing user still has full access.

statistics=/re_gecip/machine_learning/LOEUF_classifier_tool/Data_output/statistics.csv

log=/re_gecip/machine_learning/LOEUF_classifier_tool/Data_output/log.txt


chmod 700 "$statistics"

chmod 700 "$log"

printf "$COUNTER"
printf "\n"
printf "done"

# Counter will tell us the number of iteration or genomes analised.
# Handy when large number of genomes are analysed.