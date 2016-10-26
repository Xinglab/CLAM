#!/bin/bash
# -- our name ---
#$ -S /bin/bash
#$ -R y
#$ -V
# Make sure that the .e and .o file arrive in the
# working directory
#$ -cwd
#Merge the standard out and standard error to one file
#$ -j y
#
# Send mail at submission and completion of script
#$ -m be
#$ -M your@email.com
#

script_dir="./CLAM/bin"

echo "bamfile: $1"
echo "output folder: $2"
echo "tmp folder: $3"
echo "is_stranded: $4"

if [ $4 -eq 1 ]
then
echo "is stranded"
python $script_dir"/CLAM.lite_aligner.py" --verbose -o $2 -t $3 --is-stranded $1
else
echo "unstranded"
python $script_dir"/CLAM.lite_aligner.py" --verbose -o $2 -t $3 $1
fi

if [ $4 -eq 1 ]
then
python $script_dir"/CLAM.fdr_peak.MP.py" --verbose --is-stranded -o $2 -t $3 --pval-cutoff=0.001 --pval-method=2 -p CLAM_peak.bed --ThreadN=30 --max-iter=1000
else
python $script_dir"/CLAM.fdr_peak.MP.py" --verbose -o $2 -t $3 --pval-cutoff=0.001 --pval-method=2 -p CLAM_peak.bed --ThreadN=30 --max-iter=1000
fi
