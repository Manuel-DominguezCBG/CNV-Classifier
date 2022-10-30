
# CNV Classifier 
Manuel Dominguez manuel.dominguezbecerra@nhs.net
****
#### A bioinformatics tool to annotate and classify CNVs in the Genomics England (GEL) Research Environment.

##### Table of Contents 
****
* [Introduction](#introduction)  
* [Motivation](#motivation)
* [Usage](#usage)
* [Result](#results)
* [Structure](#structure)
* [Third-parties](#third-parties)
* [Troubleshooting](#troubleshooting)
* [Change control](#control)
* [Future development](#Future_development)
* [License](#license)
 
## Introduction
This software is the product of a Master's dissertation by a Clinical Bioinformatician who works in a clinical genomics lab. It aims to improve Copy Number Variation (CNV) analysis by satisfying the specific needs of a group of clinicians, clinical scientists and researchers. Since one of these users' requirements was to develop this tool in the [GEL Research Environment](https://www.genomicsengland.co.uk/research/research-environment), this application is designed to work in this High-Performance Cluster computer (Helix). Please be aware that this is not going to work on your computer unless you make the appropriate changes. This public repository is to show the main code. Some parts of the program (e.g. datasets and VCF files used for testing) are missing here because of the restrictions imposed by the GERE related to [what can be exported](https://research-help.genomicsengland.co.uk/display/GERE/What+you+can+and+can%27t+export).

### Motivation
In a healthy genome, the number of CNVs disrupting the function of Protein-Coding Gene (PCG) is large (e.g. around 100 CNVs overlapping with PCG has been seen in a 70k-genomes analysis carried out with this application). Secondly, between 70-75% of human genes have non-known functions.  Mendelian phenotypes caused by variants have been therefore identified only in approximately 3,000 PCG from a total of 20,000 PCG predicted. As a result, when disorders may involve variation in multiple genes and clinicians interrogate WGS/WES data looking for CNVs responsible for the patient phenotype, the large number of CNVs disrupting the function of genes found in the healthy population is a challenge. The aim of this project is to improve this analysis by carrying out the following tasks. First, CNVs are annotated to identify the exons, transcript and genes affected. Once we know the affected genes they are classified by the LOEUF score. This is a continuous metric that estimates the tolerance levels of genes against LoF variation. With this tool, users can prioritise genes because low LOEUF scores genes are probably haploinsufficient genes https://www.nature.com/articles/s41586-020-2308-7. 

Secondly, the application looks for compound heterozygous SNPs and indels. This is because autosomal recessive genes normally have intermediate LOEUF scores. If a CNV is affecting one copy of one gene, the identification of compound heterozygous variants is important because the second "hit" in that gene may completely “knock-out” the gene function.

CNV-Classifier is able to identify that and generates reports that allow a user (with a bioinformatic background) to carry out further analysis.  

## Usage
This is a command-line tool made mostly by Python. The main scripts are run.py and automating.sh. Both can be found at
```bash
/re_gecip/machine_learning/LOUEF_classifier_tool/Application
```


 To run one genome at a time the script to execute is run.py. I will explain later how to use automating.sh to analyse 2 or more genomes automatically. To use run.py, please follow the following steps

1. Open the terminal, access Helix HPC (CNV-classifier can be run un the loging terminal but GEL RE does not recommend to directly run applications on login nodes) and load the Bedtools module
```bash
module load bio/BEDTools/2.19.1
```

If GEL RE upgrade bedtools, this may give an error. If this occurs, check for the last available version by typing

```bash
module avail
```

2. Activate your conda environment. 
A conda environment is a directory that contains a specific set of packages previously installed. This is very important because this allows you to use the specific libraries I installed when I developed CNV-Classifier. In Helix, activate conda
```bash
source /resources/conda/miniconda3/bin/activate
```

To replicate the environment, use the .yml file created 
 ``` /re_gecip/machine_learning/LOUEF_classifier_tool/Application ```
 and run
```bash
conda env create -f requirements.yml.
```
This may take a while! After this, you need to activate the environment which name is cnv

```bash
conda activate cnv
```

If everything is ok you will see how the name of this will appear at the beginning of the line of the terminal ```(cnv_classfier)```

3. (Optional but recommended) Run the run control test.
People who maintain GEL RE may introduce changes that could affect the behaviour of this program. If this occurs (this has happened twice during the development of this program) I have no control over this. The Conda environment  partially prevents these changes but CNV-Classifier also uses modules that are outside Conda such as Bedtools or VEP. To avoid this, I have designed a Structural Variation VCF file and a Standard VCF file with manually created variants. When running this control test, this tool analyses these two files, creates the report files and verifies that the results are ok. In other words, this easy-to-use test lets you know if everything is working ok.

```bash
python run.py  run_control
```
In the terminal, if you get the message ` VALIDATION-SUCCESS: MD5 hash checking ok` that means the application should work fine. If you get an error message please contact me.

4. Run the application by selecting the path and the file name of the Structural Variation (SV) VCF file. Like this
```bash
python run.py  directory/where/genome/is/file_name.SV.vcf.gz
# This tool takes variants in a SV Variant Call format
# You don't need to decompress the file. You don't need to copy the VCF file in your personal space. 
# Genomes in GEL RE are protected (e.g. you cannot modify the original files) but this application 
# only read them and GEL allow you to do this.
```
The results are found in the folder `Data_output`. There is a section below (Results) explaining the content of the different outputs generated by this application.

5. If you wish to run many genomes. This application has been designed to be flexible and be able to generate the same results as when running only one genome. To do this, you only need to add the directory and file name of the genomes you wish to analyse in the file `'/re_gecip/machine learning/Manuel/LOEUF_tool/Resources/genome_path/For testing/automating.txt'`. One genome per line. After this, just type
```bash
bsub < automating.sh
```
run.py can be executed in the terminal of Helix with no problem. However, running many genomes exceeded the computational resources assigned to the terminal. You can execute programs in the Helix defining the resources you need using `BSUB`. Actually, you don't need to worry about configuring `BSUB` because this is done for you in `automating.sh`.  `automating.sh` is a shell script that submits a job for you satisfying all the needs that the computer needs to perform the analyziz. The same number of files are generated in the same location when you execute run.py. The results of the genomes are listed one below the other in the same order in which they are analyzed. You will find the genome identification in the first column of the file output. In conclusion, both scripts produce the same number of output files, but if you analyze >2 genomes, the results are one below the other in the same file. We have done this (as requested by the user) to make it easier to manipulate data in subsequent analyses.

`BSUB` is the command-line tool implemented in this HPC to submit a job for batch execution on hosts. You can check the stage of your submission by typing `bjob -u <your_user_name_here>`. `BSUB` generated a couple of extra files in the `Data_output`. These are `.err` and `.out`. The first informs about any errors. Sometimes, it also adds some lines with warnings or some info. If there is not an ERROR message, everything should be ok.  The .out returns value information about the resources needed, the execution time, etc. These are for advanced users, to make it easier to understand this tool, I have created the log.txt file that explains the behaviour of the tools in a friendly way.


## Results 
This application generates 4 files. 
- **` log.txt`** The historical record of usage patterns, activities and the outcome of important events are saved in this file while CNV-Classifier is running. An execution carried out correctly looks like this:

```bash
(cnv) [user@corp.gel.ac@phpgridzlogn001 Application]$ python run.py directory/where/genome/is/LP3000448-DNA_A09.SV.vcf.gz
2022-07-06 09:37:38,805 INFO:        APPLICATION STARTING
2022-07-06 09:37:39,001 INFO:        statistics file created
2022-07-06 09:37:38,806 VALIDATION:  Input is a structural VCF file
2022-07-06 09:37:38,806 INFO:        Sample running: LP0000001-DNA_A09.SV.vcf.gz
2022-07-06 09:37:39,334 VALIDATION:  File name and sample name match
2022-07-06 09:37:39,364 SUCCESS:     CNVS filtered
2022-07-06 09:37:39,364 INFO: The number of CNVS found in this VCF file are: 410
2022-07-06 09:37:39,820 SUCCESS:     CNV mapped (genes)
2022-07-06 09:37:42,207 SUCCESS:     CNV mapped (exons)
2022-07-06 09:37:42,211 VALIDATION:  Gene annotation ok? 100.0% True
2022-07-06 09:37:42,215 VALIDATION:  Exon annotation ok? 100.0% True
2022-07-06 09:37:43,959 INFO:        Looking for Compound Heterozygous Variants
2022-07-06 09:37:43,961 INFO:        Standard VCF file exists
2022-07-06 09:37:55, 295 INFO:       Starting VEP annotation (this may take a while)
2022-07-06 09:39:06,764 SUCCESS:     VEP annotation done
2022-07-06 09:39:06,794 INFO:        Number of High Confident compound heterozygous variants: 2
2022-07-06 09:39:07,132 INFO:        Application finished correctly
Doing some cleanup
```
This information is showing in the terminal if you run `run.py`. But with `run.py` and `automating.sh`, this is always recorded in the `log.txt`.  If you run automating.sh with 2 or more genomes, this information is provided for each genome. Therefore, the behaviour of the analysis is provided for each genome analysed. You can see that the name of the genome is provided so you always know what happened with what particular analysis. If errors happen, this will be also found in this file and you will know in what part of the workflow the error occurs. I created this because when using `BSUB` there is a .err file however this provided limited info about Python errors.

- **` The CNV_annotated_by_genes.tsv`**
This is the main report generated by this bioinformatic application. This report shows the affected genes one per line. Let's see this table with a couple of genes as examples.

| Participant identification | Gene                   | CNV_type_&_QUAL | Canonical | MANE | LOEUF metrics (max) | LOEUF metric (min) | HGVSc                            | GT      | Consequence | Uploaded_variation | Allele | Phenotype description     | CNV overlapping with gene |
|---------------------------|------------------------|-----------------|-----------|------|---------------------|--------------------|----------------------------------|---------|-------------|--------------------|--------|---------------------------|---------------------------|
| LP00000001-DNA_A01        | MAP2K3-ENSG00000034152 | GAIN-21         | Yes       | Yes  | 1.636               | 0.65               | HGVSc=ENST00000613338.4:c.217C>T | ZYG=HET | stop_gained | rs55796947         | T      |                           | 25883.0                   |
| LP00000001-DNA_A01        | DGCR6-ENSG00000183628  | LOSS-23         | Yes       | Yes  | 1.936               | 1.619              |                                  |         |             |                    |        | VELOCARDIOFACIAL SYNDROME | 6060.0                    |
(This and the following data is always **fake data**)

Table description by columns

-- **Participant identification** This is how GEL identify the participants

-- **Gene**: Gene name and Ensembl Gene ID (Ensembl gene 107, GRCh38.p13)

-- **CNV_type_&_QUAL**: This information is taken from the SV VCF file. Structural variation is taken from the SV VCF file. This file record a large amount (~15.000) of structural variations. Two variant callers (CANVA and MANTAS) are used by the GERE to call these variants. In this project, we are only interested in a small fraction of this structural variation. Users are interested on TYPE=CNV from the CANVA algorithm (type only can be LOSS or GAIN). As mentioned above, this is the main reason this application look for satisfying the very specific needs of the members of my team and likely you may need to introduce minor changes to adapt this tool. This column also shown the quality of the variant that is affecting the gene. This is the Phred-scaled quality score for the assertion made [(The Variant Call Format Version 4.2 Specification)](https://samtools.github.io/hts-specs/VCFv4.2.pdf).

-- **Canonical**: Yes if the CNV affects the canonical transcript of that gene, otherwise NO. Canonical transcripts according to Ensembl gene 107, GRCh38.p13.

-- **MANE**: Same criteria for MANE.  It should be noticed that this application only works with GRCh38 as 95% of the Genomes analysed in GEL are in this Genome Reference. If you want to make this work with the previous version you may need to modify the datasets.

-- **LOEUF metrics (max)**. LOEUF is a transcript-specific metric. So, a gene may have more than one transcript. Normally, canonical ones present a lower value. This and the following column show you the highest and the lowest metrics associated with that gene. The transcripts affected are provided in the following report explained below.

-- **LOEUF metrics (min)**

-- **HGVSc**: One of the main requirements in this project was to identify compound heterozygous variants (short variants) in genes affected by CNVs. We have seen not many compounds with heterozygous variants (~2)  per genome. How these variants are filtered is explained below. At the moment, variants found in this column are only High Confident (HC) predicted Loss-of-functions (LoF) SNPs and Indels found in the exons of transcripts affected by a CNV. This column tells you the HGVSc nomenclature.

-- **GT**: The Genotype of this variant. HOM or HET

-- **Consequence** Three possible consequences based on the filtering process: Frameshift_variant, Stop_gained_variants or Splice_donor_variant. This is because we have used [LOFTEE](https://github.com/konradjk/loftee) to filter and select HC predicted LoF variants. 

-- **Uploaded_variation**: If this variant is found in the dbSNP database (VEP 99), this is reported in this column. I recommend checking this because new variants are added daily and working in GEL I cannot connect this with the most recent version.

-- **Allele**. For SNPs, this tells you the ALT. Not really needed as this is given by HGVSc. So, this may be deleted in future versions.

-- **Phenotype_description**: Phenotypes such as diseases and traits are associated with the genes in the report. This information is taken from  [COSMIC](http://www.sanger.ac.uk/genetics/CGP/cosmic/), [ClinVar](http://www.ncbi.nlm.nih.gov/clinvar/), [dbGaP](http://www.ncbi.nlm.nih.gov/gap), [GIANT](http://www.broadinstitute.org/collaboration/giant/index.php/Main_Page), [HGMD-Public ](http://www.magicinvestigators.org/), [MAGIC](http://www.sanger.ac.uk/genetics/CGP/cosmic/) and [NHGRI-EBI GWAS catalog](http://www.ebi.ac.uk/gwas/).

-- **CNV_overlapping_with_gene**. Clinicians requiered to know the number of bases in the gene affected by the CNV. It should be more informative if the proportion is given instead. (For future development).


- **`The CNV_annotated_by_transcripts`**
This is the complementary report and provides the full data about the CNV analysis carried out by this program. This report is very large because this annotated CNV instead genes. For example, if a CNV affects 8 exons of a gene involving 2 transcripts, 8 rows will be provided showing information about these exons, the transcripts and so on.

| Header                   |Description                    | Example |
|-------------------------|--------------------------------------------------------------| --- |
| CHROM                   | Chromosome where CNV identified (taken from the VCF file) [String]  | Chr1
| START                   | Start coordinate of CNV (taken from CNV file) [Interger]            | 152582182 |
| END                     | End coordinate of CNV (taken from CNV file) [Interger]              | 152615293 |
| CNV_TYPE                | Delition or duplication + QUAL (taken from VCF file)    [String]    | LOSS-15
| CNV_ID                  | ID provided by the VCF. This is participant ID + coordinates        |LP0000001-DNA_A1_G12|CANVAS:LOSS|...
| Chromosome_g            | Chromosome of the gene (from gene annotation dataset)  [String]      | Chr1
| Transcript_start_g      | Start coordinates of the transcripts      [Interger]                 | 152613811
| Transcript_end_g        | End Coordinates of the transcripts      [Interger]                   | 152613811
| Transcript_ID           | Ensembl transcript ID  [String]                                      | ENST00000335633, if not gene affected = NO_GENE_AFFECTED
| Canonical_g             | Chromosome of the gene taken from genes dataset   [Boolean]           | TRUE, if not gene affected = -1
| Gene_stable_ID          | Ensembl Gene ID                [String]                              | ENSG0000018723, if not gene affected = NO_GENE_AFFECTED
| Gene_name               | Gene name                    [String]                                | LCE3B
| amount_overlap          | Number of bases that overlap between the gene and the CNV [Interger]    | 1234567879, if not gene affected = 0
| Chromosome_e            | Chromosome of the affected exon (from gene annotation dataset)  [String]| Chr1
| Exon_start              | Start coordinate of exon (from gene annotation dataset)  [Interger]            | 152613811 if not gene affected = -1
| Exon_end                | End coordinate of exon (from gene annotation dataset)     [Interger]            | 152614098 if not gene affected = -1
| Exon_rank_in_transcript | Exon rank (from gene annotation dataset)  [Interger]                             | 1
| Exon_stable_ID          | Ensembl exon ID (from gene annotation dataset)        [String]               | ENSE0000133845 if not exon affected = -1
| Constitute_exon_e       | 0 = no, 1 = yes (from gene annotation dataset)     [Interger]                 | 1
| Transcript_stable_ID    | Ensembl transcript ID taken from exon dataset                | ENST00000335633
| Gene_name_e             | Gene name taken (from exon annotation dataset)                          | LCE3B
| oe_lof_upper            | LOEUF metrics (continuous 0 to 2)          [String]                  |  1.425, if not gene affected = .
| oe_lof_upper_bin        | LOEUF metrics is also given by bins 0 to 9                   | 9, if not gene affected = .
| pass_validation_exon    | TRUE/FALSE   [Boolean]                                       | TRUE



- **` The statistic.csv file`**

This file provides a dataset with information related to the CNVs analysed. This stems from a suggestion from one of the researchers on our team. It was proposed to create a function that would allow the collection of data related to genomes, SV, and affected genes. During the development of this program, I have incorporated this function to collect information that I personally considered important to later perform a () analysis. The idea was to create a simple function that would allow the user to modify it according to their own interest. As it is developed CNV-Classifier, the table created  looks like followed

| File name | Nof variants in SV VCF file | Nof of CNVs | Nof Loss | Nof Gain | Nof CNVs no affecting CPG | Nof CNVs no affecting CPG (Gain) | Nof CNVs no affecting CPG (Loss) | Nof transcripts affected | Nof CNVs affecting 1 G | Nof CNVs affecting 2 G | Nof CNVs affecting 3 G | Nof CNVs affecting 4 G | Nof CNVs affecting >5 G | How likely CNVs affect >1 transcripts (%) | Nof G affected by CNVs |
|-----------|-----------------------------|-------------|----------|----------|---------------------------|----------------------------------|----------------------------------|--------------------------|------------------------|------------------------|------------------------|------------------------|-------------------------|-------------------------------------------|------------------------|
| LP3000001 | 14471                       | 34          | 19       | 15       | 25                        | 12                               | 13                               | 42                       | 10                     | 5                      | 2                      | 0                      | 0                       | 41.5                                      | 17                     |
| LP3000001 | 14531                       | 37          | 20       | 17       | 30                        | 15                               | 15                               | 40                       | 11                     | 8                      | 1                      | 1                      | 0                       | 45.5                                      | 20                     |

This function creates this file if there is not one found and put the header. Then, for each new genome analysed, the function adds the data to one line for each genome analysed. The result is a dataset ready to be used by the researcher.

## Structure 
```
├── Application
│   ├── automating.sh   
│   ├── classes
│   │   ├── SNP_annotation.sh
│   │   ├── __init__.py
│   │   ├── add_metrics.py
│   │   ├── annotation.py
│   │   ├── automated_testing.py
│   │   ├── control_sample.py
│   │   ├── data_presentation.py
│   │   ├── data_presentation2.py
│   │   ├── filter.py
│   │   ├── import_data.py
│   │   ├── ini_config.py
│   │   ├── link.py
│   │   ├── log.py
│   │   ├── main.py
│   │   ├── read_file.py
│   │   └── short_variants.py
│   ├── config.ini
│   ├── requiriments.yaml
│   └── run.py
├── Data_output
│   ├── RESULTS
│   └── tmp
├── README.md
├── test
│   └──  test_application.py
├── Resources   
│   ├── For_annotation
│   └── genome_path
│       ├── ALL_GRCh38
│       ├── All_RareDiseaseGermline_GRCh38
│       └── For_testing
└── run_control          
    ├── MD5_checking
    └── manually_created_VCF
```
- The datasets needed to run this program are not in this repository due to the size of these files.
- The run control analysis is not provided in this repository because they have fake but also real variants from GEL participants and we cannot export this outside GEL RE.
- The config.ini file allows you to remove variants based on a cutoff. By default, CNVs with QUAL lower than 14 are not considered in the final report. By default, only CNVs that PASS all filters are taken. 


## Third-parties

- Bedtools (v2.19.1)
- BCFtools (v1.9)
- VEP (v99.1)
- LOFTEE (v1.0.2)

## Troubleshooting 
The main issues you may find when working with this application are related to the tools this program uses.

1. ```Error: Module not found``. This application uses Bedtools and VEP. For a reason it is unclear, sometimes the HPC takes more time than normal to load these modules. Based on GEL documentation, before using VEP, you have to unload any other module you have loaded. So, when analysing many genomes this causes huge delays. Sometimes you would need to cancel (to do this `bkill <your_user_name>` this will cancel all jobs you have submitted) and resubmit. Or sometimes the module is not found. If not found, make sure its name has not been changed due to a new version. For reasons like this, running the test will let you know if everything is ok.

2. ```Error: Unable to open file ... Exiting``` This happens with Bedtools when parsing the standard VCF (the one with short variants)  the error message is not very informative and that not show a reason. However, I have solved it by splitting the VCF file which is telling me that the error is due to the file being too large to be analysed at once. After running all genomes found in the last release of GEL, I have found that this happens very often. So, the application is able to detect this error and takes an alternative to generate the final report. Genes will be annotated with CNVs but no compound heterozygous variants are added. To detect this error, you have two alternatives. This is notified in the log file or in the statistic file, in the column that informs the number of compounds heterozygous, you will find `0 due to an error`. For the next development, if this error is found, I would like to implement something that split the VCF file into two files, then analyse both, merge data and provide the compound heterozygous variants. 
3. Recently, the GERE is doing maintenance work and sometimes CONDA does not work. We have also noticed that although sometimes CONDA works, it does not allow you to use env. previously created by us or does not let you create a new env. As a preventive measure (if this happens again), the "SigProfilerExtractor-py3-test" env. has all the necessary libraries to make CNV-Classifier work correctly. Tested this at 20/Oct/22. 

## Change control

### CNV-Classifier v1.1.0-a.1
The program is able to do CNV annotation (genes and transcripts).
### CNV-Classifier v1.1.0-a.2
The program is able to compound heterozygous annotation. Log file and statistics dataset features added. Config file allow now the user to modify QUAL and FILTER features.
### CNV-Classifier v1.1.0-a.3
Main and complementary reports are now created as the user requested. Automating.sh allow the user to analyse >1 genome.
### CNV-Classifier v1.1.0-rc.3
Minor changes before Release. Added The control-run test.
### CNV-Classifier v1.1.0
First release. 

## Future development 

Although there is a .config file that allows selection of CNVs based on quality, this program is not flexible enough to allow selection of different types of variants (e.g. complex variants). It would be interesting to allow the user to select the type or variant types via command line options (e.g. Python run.py --CNV_only, --Del_only and so on) or by using the .config file. **TODO** The selection of the variants has to be preferably at the beginning of the program (as it currently happens) because if we carry out the entire analysis with the entire ST VCF file (~15k variants), the program slows down a significatly. Better filtering them at the beginning.

## License
Copyright 2022 <Manuel Dominguez>

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

