# CLAM Version 1.1.3
# CLIP-seq Analysis of Multi-mapped reads

## Table of Contents
 - [Introduction](#introduction) 
 - [Installation](#installation)
 - [Input](#input)
 - [Usage](#usage)
   - [preprocessor](#clam-preprocessor)
   - [realigner](#clam-realigner)
   - [peakcaller](#clam-peakcaller)
   - [permutation_callpeak](#clam-permutation_callpeak)
 - [Output](#output)
 - [Testing data](#testing-data)
 - [Contacts](#contacts)

## Introduction
CLAM is a general toolkit for re-aligning multi-mapped reads in CLIP/RIP-seq data and calling peaks.

For details, please read our [NAR paper](https://academic.oup.com/nar/article/45/16/9260/4077049/CLIP-seq-analysis-of-multi-mapped-reads-discovers).

[TOC](#clip-seq-analysis-of-multi-mapped-reads)

## Installation
CLAM v1.1 works under Python 2. Please click and download the latest version from the releases. Once unzip the file, type
```
$ python setup.py install
```
in your terminal and this will automatically install CLAM in your currently working python.

You should have already installed "pysam" using pip/conda for your python interpreter.
If not, you can check the detailed requirements in the file "requirements.txt", or type
```
pip install -r requirements.txt
```
to install those requirements manually.

[TOC](#clip-seq-analysis-of-multi-mapped-reads)

## Input
The input for CLAM is a sorted or unsorted BAM file of CLIP-seq alignment and a gene annotation file in GTF format.

In the case of RIP-seq or eCLIP, a BAM file for IP experiment and a BAM file for 
Control/input experiment are taken together as input.

![#f03c15](https://placehold.it/15/f03c15/000000?text=+) Note: As in the released version v1.1.0, the `read_gtf` function had a bug and required the Gencode format GTF (i.e. last column of GTF
matches gene_id "(xx)" ) to proceed the peak calling. This bug has been fixed in the github repository and has been fixed in later releases/patches.

[TOC](#clip-seq-analysis-of-multi-mapped-reads)

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

Below we briefly describe what each subcomand does and provide an example command-line run.

#### CLAM preprocessor
This subcommand (new v1.1) will prepare the input files for CLAM pipeline. As of the current version (v1.1), it looks for 
reads passing QC, splits the input bam file by sorting them into `unique.sorted.bam` and `multi.sorted.bam`, 
and adding an additional tag "RT" (short for Read Tag) to each alignment based which read tagger function the user supplied.

Note that you can also run `CLAM realigner` directly, which will call `preprocessor` and automatically determine
if `preprocessor` has been called in the output folder. 

If you don't want to run `realigner`, you can also run `peakcaller` directly after `preprocessor`.

Example run:
```
CLAM preprocessor -i path/to/input/Aligned.out.bam -o path/to/clam/outdir/ --read-tagger-method median
```

#### CLAM realigner
This subcommand will run expectation-maxmization to assign the multi-mapped reads in a probablistic framework. 
More details about the EM model is described in our NAR paper.

Note when `--retag` is specified, `realigner` will re-run `preprocessor` regardless; otherwise, it will use 
the prepared files in `outdir` if available.

Example run:
```
CLAM realigner -i path/to/input/Aligned.out.bam -o path/to/clam/outdir/ --read-tagger-method start --retag
```

#### CLAM peakcaller
This subcommand (new in v1.1) will call peaks by looking for bins enriched with IP reads over control, specifying a 
Negative-binomial model on observed read counts.

Note we can specify both `unique.sorted.bam` (from `preprocessor`) and `realigned.sorted.bam` (from `realigner`) and 
separte the two file paths by a space, to call peaks using the combination of uniquely- and multi-mapped reads.

Alternatively, we can also only input `unique.sorted.bam`; this will allow CLAM to call peaks using only uniquely-
mapped reads.

Example run:
```
CLAM peakcaller -i path/to/IP/outdir/unique.sorted.bam path/to/IP/outdir/realigned.sorted.bam \
-c path/to/CTRL/unique.sorted.bam path/to/CTRL/realigned.sorted.bam \
-o path/to/peaks/outdir --unstranded --binsize 100 \
--gtf path/to/gencode.v19.annotation.gtf
```

#### CLAM permutation_callpeak
This subcommand will call peaks using permutation by randomly placing reads along the gene.
More details about the permutation procedure is described in our NAR paper.

Example run:
```
CLAM permutation_callpeak -i path/to/outdir/unique.sorted.bam path/to/outdir/realigned.sorted.bam \
-o path/to/peaks/outdir -p 8 \
--gtf path/to/gencode.v19.annotation.gtf
```

[TOC](#clip-seq-analysis-of-multi-mapped-reads)

## Output
The output of the re-aligner is "realigned.sorted.bam" (previously "assigned_multimapped_reads.bam" in v1.0), 
which is a customized BAM file following SAM format. 
Note that the re-aligned weights are stored in "AS:" tag, so please be aware and do not change/omit it.
Output of re-aligner could also be seen as an intermediate file for CLAM pipeline.

The output of the peak-caller is a bed file following NarrowPeak format. It is a 10-column [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) format file. 

If you run permutation peak caller (as in v1.0), there will be only one output file called "narrow_peak.permutation.bed".
Hence a peak with "combined" but no "unique" on the fifth column indicates this is a rescued peak; both "unique" and 
"combined" as common peak; or lost peak otherwise.

If you run model-based peak caller (new in v1.1), depending on the specified paramter (whether you turned on `--unique-only`), 
the output will be either "narrow_peak.unique.bed" for peaks called using only uniquely-mapped reads; or 
"narrow_peak.combined.bed" for peaks called when adding realigned multi-mapped reads.

[TOC](#clip-seq-analysis-of-multi-mapped-reads)

## Testing data
Once downloaded the CLAM source code, please download the hnRNPC iCLIP dataset from [here](https://xinglab.cass.idre.ucla.edu/public/zijun/CLAM/test_data/hnRNPC_iCLIP_rep1_E-MAT-1371_novoalign.sorted.bam).

Then run CLAM by using `realigner` and `permutation_peakcaller` on the dataset; if finished correctly, you should have rescued peaks at these two loci:

chr11:82,624,179-82,626,008

chr20:37,054,180-37,055,310

[TOC](#clip-seq-analysis-of-multi-mapped-reads)


## Contacts
Zijun Zhang [zj.z@ucla.edu](mailto:zj.z@ucla.edu)

Yi Xing [yxing@ucla.edu](mailto:yxing@ucla.edu)

If you found a bug or mistake in this project, we would like to know about it. 


[TOC](#clip-seq-analysis-of-multi-mapped-reads)

## Copyright and License Information
Copyright (C) 2016-2017 University of California, Los Angeles (UCLA) Zijun Zhang and Yi Xing

Authors: Zijun Zhang and Yi Xing

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see [http://www.gnu.org/licenses/](http://www.gnu.org/licenses/).
