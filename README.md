# CLAM Version 1.1.0
# CLIP-seq Analysis of Multi-mapped reads

## Download the latest version [here](https://github.com/Xinglab/CLAM/releases/download/1.1.0/CLAM_v1.1.zip).

## Introduction
CLAM is a general toolkit for re-aligning multi-mapped reads in CLIP/RIP-seq data and calling peaks.

## Installation
CLAM v1.1 works under Python 2. Please click and download the latest version from the releases. Once unzip the file, type
```
$ python setup.py install
```
in your terminal and this will automatically install CLAM in your currently working python.

You should have already installed "pysam" using pip/conda for your python interpreter.
If not, you can check the detailed requirements in the file "requirements.txt", or type
```
pip -r requirements.txt
```
to install those requirements manually.

## Usage
CLAM is run through issueing subcommands. Currently there are four subcommands available:
preprocessor, realigner, peakcaller and permutation_callpeak.
```
usage: CLAM [-h] [--version]
            {preprocessor,realigner,peakcaller,permutation_callpeak} ...

CLAM -- CLip-seq Analysis of Multi-mapped reads

positional arguments:
  {preprocessor,realigner,peakcaller,permutation_callpeak}
    preprocessor        CLAM Preprocessor: tag read alignments to specific
                        locus
    realigner           CLAM Realigner: realign multi-mapped reads using
                        expectation-maximization
    peakcaller          CLAM Peakcaller: negative binomial model-based peak
                        calling combining unique- and multi-reads
    permutation_callpeak
                        CLAM permutation peakcaller: call peaks using
                        permutation (as in v1.0.0)

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit

For command line options of each sub-command, type: CLAM COMMAND -h
```


## Output
The output of the re-aligner is "assigned_multimapped_reads.bam", which is a customized BAM file following SAM format. Note that the re-aligned weights are stored in "AS:f" tag, so please be aware and do not change/omit it.
Output of re-aligner could also be seen as an intermediate file for CLAM pipeline.

The output of the peak-caller is a bed file following NarrorPeak format. It is a 10-column [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) format file, separated by tabulate and ordered as 
Hence a peak with "combined" but no "unique" on the fifth column indicates this is a rescued peak; both "unique" and "combined" as common peak; or lost peak otherwise.

## Testing data
Once downloaded the CLAM source code, please download the hnRNPC iCLIP dataset from [here](http://www.mimg.ucla.edu/faculty/xing/CLAM/hnRNPC_iCLIP_rep1_E-MAT-1371_novoalign.sorted.bam).

Then run CLAM on the dataset; if finished correctly, you should have rescued peaks at these two loci:

chr11:82,624,179-82,626,008

chr20:37,054,180-37,055,310



## Contacts
Yi Xing [yxing@ucla.edu](mailto:yxing@ucla.edu)

Zijun Zhang [zj.z@ucla.edu](mailto:zj.z@ucla.edu)

If you found a bug or mistake in this project, we would like to know about it. Before you send us the bug report though, please check the following:

1. Are you using the latest version? The bug you found may already have been fixed.
2. Check that your input is in the correct format and you have selected the correct options.
3. Please reduce your input to the smallest possible size that still produces the bug; we will need your input data to reproduce the problem, and the smaller you can make it, the easier it will be.

## Copyright and License Information
Copyright (C) 2016-2017 University of California, Los Angeles (UCLA) Zijun Zhang and Yi Xing

Authors: Zijun Zhang and Yi Xing

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see [http://www.gnu.org/licenses/](http://www.gnu.org/licenses/).
