# ToxinPred3.0
A method for predicting toxicity of the peptides
# Introduction
ToxinPred3.0 is developed for predicting, mapping and scanning toxic/non-toxic peptides. It uses only composition based features for predicting toxic/non-toxic peptides. The final model also deploys a motif-based module which has been implemented using MERCI. More information on ToxinPred3.0 is available from its web server http://webs.iiitd.edu.in/raghava/toxinpred3. Please read/cite the content about toxinpred3.0 for complete information including algorithm behind the approach.

## Create an Environment
1- Install miniconda. You need it to create an environment and install the packages this project need. for more information, visit https://www.anaconda.com/docs/getting-started/miniconda/main

use the following command in your terminal to create an environment with python.

```
conda create -n toxinpred-env python=3.8
```

After installation was complete, use the following command to active (go to) the environment.

```
conda activate toxinpred-env
```

Now, your terminal knows python3, version 3.8 and now you can run the program, but first, lets install some packages.

## Installing packages

In order to install all the packages you need, you can run the following command, or install the packages one by one, if you know how, remember that the project is version sensitive. Hence the right package with wrong version, will not work. The following command, install the right packages with the right version, so you do not have to worry about anything.

```
pip install -r requirements.txt
```

## Minimum USAGE
To know about the available option for the standalone, type the following command in your terminal after all the packages are installed. remember that the environment you made with conda should be activated, otherwise, you will hit a not found name error:
```
toxinpred3.py -h
```
To run the example, type the following command:
```
toxinpred3.py -i peptide.fa

```
**Full Usage**: 
```
Following is complete list of all options, you may get these options
usage: toxinpred3.py [-h] 
                     [-i INPUT]
                     [-o OUTPUT]
                     [-t THRESHOLD]
                     [-m {1,2}] 
                     [-d {1,2}]
```
```
Please provide following arguments

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input: protein or peptide sequence in FASTA format or
                        single sequence per line in single letter code
  -o OUTPUT, --output OUTPUT
                        Output: File for saving results by default outfile.csv
  -t THRESHOLD, --threshold THRESHOLD
                        Threshold: Value between 0 to 1 by default 0.38
  -m {1,2}, -- model Model
                        Model: 1: ML model, 2: Hybrid model, by default 2
  -d {1,2}, --display {1,2}
                        Display: 1:Toxin peptide, 2: All peptides, by
                        default 1

```

**Input File**: It allow users to provide input in two format; i) FASTA format (standard) (e.g. peptide.fa) and ii) Simple Format. In case of simple format, file should have one peptide sequence in a single line in single letter code (eg. peptide.seq). 

**Output File**: Program will save result in CSV format, in case user do not provide output file name, it will be stored in outfile.csv.

**Threshold**: User should provide threshold between 0 and 1, please note score is proportional to toxic potential of peptide.

**Models**:  In this program, two models have been incorporated;  i) Model1 for predicting given input peptide sequence as toxic and non-toxic peptide using Extra tree based on amino-acid composition (AAC) and di peptide composition (DPC) of the peptide; 

ii) Model2 for predicting given input peptide sequence as toxic and non-toxic peptide using Hybrid approach, which is the ensemble of Extra tree + MERCI. It combines the scores generated from machine learning (ET), and MERCI as Hybrid Score, and the prediction is based on Hybrid Score.


ToxinPred3.0 Package Files
=======================
It contain following files, brief description of these files given below

INSTALLATION  	: Installation instructions

LICENSE       	: License information

merci : This folder contains the program to run MERCI

README.md     	: This file provide information about this package

toxinpred3.py 	: Main python program

peptide.fa	: Example file contain peptide sequences in FASTA format

peptide.seq	: Example file contain peptide sequences in simple format

## Installation via PIP
User can install ToxinPred3 via PIP also
```
pip install toxinpred3
```
## Reference: 
Rathore AS, Arora A, Choudhury S, Tijare P, Raghava GPS (2024) ToxinPred3.0:An improved method for predicting the toxicity of peptides. 
Comput Biol Med. 179:108926 . https://doi.org/10.1016/j.compbiomed.2024.108926

Rathore AS, Arora A, Choudhury S, Tijare P, Raghava GPS. ToxinPred3.0:An improved method for predicting the toxicity of peptides. bioRxiv 2023.08.11.552911; doi: https://doi.org/10.1101/2023.08.11.552911
