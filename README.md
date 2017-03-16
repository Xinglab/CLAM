# CLAM Version 1.0.0
# CLIP-seq Analysis of Multi-mapped reads

## Requirements
Download the latest version [here](https://github.com/Xinglab/CLAM/releases/download/1.0.0/CLAM-v1.zip).

CLAM is a two-stage algorithm implemented in Python. It is intended to be used in Unix-based environment. It has been tested on Linux with Python 2.7.3.

CLAM depends on several commonly-used Python libraries, including [pysam](http://pysam.readthedocs.io/en/latest/) and [pybedtools](https://daler.github.io/pybedtools/index.html).

A detailed dependency requirements (with version info) can be found in "requirements.txt". Alternatively, just run
```
pip install -r requirements.txt
```
and you will be good to go.

## Usage
We provide a general shell script wrapper that runs the whole pipeline sequentially with default parameters for CLIP-seq. You only need to give the paths to input bam file and output folder, and a binary (0/1) indicator for strandness:
```
$ sh runCLAM_git.sh $bam $output_dir $temp_dir $is_stranded
```
..and the CLAM pipeline's output will be generated in $output_dir as specified. Check "Output" section below to understand the file formats in the CLAM output folder.


Alternatively, if you want to dig more into the parameters, you can run the pipeline with "--help" in command line and check the options. The following should be printed to the screen:

For CLAM re-aligner,
```
$ python CLAM.lite_aligner.py --help
Usage: CLAM.lite_aligner.py <options> input_file.bam

Options:
  -h, --help            show this help message and exit
  -o OUTPUT_DIR         Output file folder [Default ./out_CLAM]
  -t TMP_DIR            Temporary file folder [Default ./tmp_CLAM]
  -w WINDOW_SIZE        Local window size for EM [Default: 50]
  --max-multihits=MAX_MULTIHITS
                        Discard reads mapped to more than <max_multihits>
                        locations. [Default: 100]
  --min-unique-reads=MIN_UNIQUE_READS
                        Discard genomic regions with less than
                        <min_unique_reads> of unique reads. [Default: 0]
  --is-stranded         Indicates if the reads are mapped with strand
                        information. [Default: False]
  --resume              Resume mode - skipping pre-processing [Default: False]
  --verbose             Verbose mode - print out all intermediate steps
                        [Default: False]
  --max-gap=MAX_GAPS    Maximum distance allowed in grouping reads. [Default:
                        50]
```

For CLAM peak-caller,
```
$ python CLAM.fdr_peak.MP.py --help
Usage: CLAM.fdr_peak.MP.py <options>

Options:
  -h, --help            show this help message and exit
  --resume              Resume mode - skipping pre-processing [Default: False]
  --verbose             Verbose mode - print out all intermediate steps
                        [Default: False]
  -o OUTPUT_DIR         Output file folder [Default ./out_CLAM]
  -t TMP_DIR            Temporary file folder [Default ./tmp_CLAM]
  -p PEAK_FILE          Output peak calling filename; if None then do not call
                        peaks  [Default none]
  --is-stranded         Indicates if the reads are mapped with strand
                        information. [Default: False]
  --extend=EXTEND       Extend to given nucleotides symmetrically at peak
                        calling [Default: 50]
  --pval-cutoff=PVAL_CUTOFF
                        Corrected p-value threshold at peak calling [Default:
                        0.001]
  --merge-size=MERGE_SIZE
                        merging window size at peak calling [Default: 50]
  --max-iter=MAX_ITER   maximum iterations for permutation tests [Default:
                        1000]
  -g GTF                GTF file [Default: ./GTF/hg19_ensembl.sorted_gene.bed]
  --ThreadN=NB_PROC     Number of threads when doing permutations. [Default:
                        4]
  --seed=SEED           Random seed for permutations. [Default: 100]
  --merge-method=MERGE_METHOD
                        Peak merging method. 1: Narrow peak 2: Broad peak
                        [Default: 1]
  --pval-method=CORRECTION_METHOD
                        Multiple testing correction method. 1: Bonferroni 2:
                        BH FDR [Default: 1]
  --call-transcriptome  Call peaks on transcriptome instead of genes with
                        multi-mappers. [Default: False]
```
And you can specify your own parameters accordingly. For example, for **m6A RIP-seq**, window size parameter (-w) and maximum gaps (--max-gap) for re-aligner should be set to 100.

## Output
The output of the re-aligner is "assigned_multimapped_reads.bam", which is a customized BAM file following SAM format. Note that the re-aligned weights are stored in "AS:f" tag, so please be aware and do not change/omit it.
Output of re-aligner could also be seen as an intermediate file for CLAM pipeline.

The output of the peak-caller is a bed file whose name is specified by user. It is a 6-column [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) format file, separated by tabulate and ordered as 
```
chr    start    end    height;fdr;gene    unique/combined    strand
```
Hence a peak with "combined" but no "unique" on the fifth column indicates this is a rescued peak; both "unique" and "combined" as common peak; or lost peak otherwise.

## Testing data
Once downloaded the CLAM source code, please download the hnRNPC iCLIP dataset from [here]()
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
Copyright (C) 2016 University of California, Los Angeles (UCLA) Zijun Zhang and Yi Xing

Authors: Zijun Zhang and Yi Xing

This program is licensed with commercial restriction use license. Please see the attached LICENSE file for details.
