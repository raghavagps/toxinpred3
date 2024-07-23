# ToxinPred3.0
A method for predicting toxicity of the peptides
# Introduction
ToxinPred3.0 is developed for predicting, mapping and scanning toxic/non-toxic peptides. It uses only composition based features for predicting toxic/non-toxic peptides. The final model also deploys a motif-based module which has been implemented using MERCI. More information on ToxinPred3.0 is available from its web server http://webs.iiitd.edu.in/raghava/toxinpred3. Please read/cite the content about toxinpred3.0 for complete information including algorithm behind the approach.

## PIP Installation
PIP version is also available for easy installation and usage of this tool. The following command is required to install the package 
```
pip install toxinpred3
```
To know about the available option for the pip package, type the following command:
```
toxinpred3 -h
```

# Standalone

Standalone version of ToxinPred3.0 is written in python3 and the following libraries are necessary for a successful run:

- scikit-learn
```
 !pip install scikit-learn==1.0.2
```
- Pandas
- Numpy

# Important Note

- Due to large size of the model file, we have compressed model. 
- It is crucial to unzip the file before attempting to use the code or model. The compressed file must be extracted to its original form for the code to function properly.



**Minimum USAGE** 

To know about the available option for the standalone, type the following command:
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

pip install toxinpred3

## Reference: 
Rathore AS, Arora A, Choudhury S, Tijare P, Raghava GPS (2024) ToxinPred3.0:An improved method for predicting the toxicity of peptides. 
Comput Biol Med. 179:108926 . https://doi.org/10.1016/j.compbiomed.2024.108926

Rathore AS, Arora A, Choudhury S, Tijare P, Raghava GPS. ToxinPred3.0:An improved method for predicting the toxicity of peptides. bioRxiv 2023.08.11.552911; doi: https://doi.org/10.1101/2023.08.11.552911
