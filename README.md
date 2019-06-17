# DASE-AG

conditional-specific Differential Alternative Splicing Events estimation method for Around-Gap regions

## Description

DASE-AG is a method for comprehensive detection of conditional-specific alternative splicing events based on de novo transcriptome assembly.

## Features

DASE-AG does not require reference genome. Instead, DASE-AG requires read files produced through RNA-seq experiments. 

## Requirement 

Python (>= 2.7)
MAFFT
NumPy
Trinity (>= v2.7)
Bowtie2
RSEM or Kallisto

## Preparation

Before execution of DASE-AG, a configuration file and a file for describing all the locations of read files are required. Check "config_example.txt" and "samples_example.txt" and adjust to your environment.

## Usage

    python dase.py -cf config.txt -out output.dat

## Author

Kouki Yonezawa (k_yonezawa@nagahama-i-bio.ac.jp)

## Reference of DASE-AG

TBA
