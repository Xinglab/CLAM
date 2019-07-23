import sys
import os
import pybedtools
import argparse as ap
import logging
from . import download_data, config

'''
Assign peaks to genomic regions
Zijun Zhang
8.1.2018
10.25.2018: wrapped to a function with document

DWK
modified to output annotation file
6.12.2019
'''

# pylint: disable-msg=too-many-function-args
# pylint: disable-msg=unexpected-keyword-arg


def parser(args):
    """DOCSTRING
    Args
    Returns
    """
    try:
        peak_in = args.peak_in
        genome = args.genome
        out_file = args.out_file
        if 'CLAM_DAT' not in os.environ or not download_data.check_genome_data(genome):
            print("Unable to locate CLAM data folder for genomic regions, will try to download.")
            print("Downloading...")
            download_data.download_genome(genome)
        genome_data = os.environ['CLAM_DAT']
        intersect_gtf_regions(
            peak_in, out_file, os.path.join(genome_data, genome))
    except KeyboardInterrupt():
        sys.exit(0)


def intersect_gtf_regions(peak_fp, outfn, gtf_dir):
    '''function: intersect_gtf_regions(peak_fp, outfn, gtf_dir)
    Intersect a peak BED file with a list of genomic region annotations (e.g. start/stop codon, UTR, intron),
    output the peak-region annotations.
    :param peak_fp: filepath to a BED-format peakquit
    :param outfn: filepath to output count file, has to end with ".txt"; annotation will be "NNN.annot.txt"

    '''
    # input arguments

    # make pybedtools objects
    print("Loading peaks...")
    peaks = pybedtools.BedTool(peak_fp)
    print("Peak file loaded.")
    print("Loading genome annotation...")
    ref_dict = {
        'exon': pybedtools.BedTool(os.path.join(gtf_dir, 'exons.bed')),
        '3UTR': pybedtools.BedTool(os.path.join(gtf_dir, '3UTRs.bed')),
        '5UTR': pybedtools.BedTool(os.path.join(gtf_dir, '5UTRs.bed')),
        'cds': pybedtools.BedTool(os.path.join(gtf_dir, 'cds.bed')),
        'intron': pybedtools.BedTool(os.path.join(gtf_dir, 'introns.bed')),
        'proximal200': pybedtools.BedTool(os.path.join(gtf_dir, 'proximal200_intron.bed')),
        'proximal500': pybedtools.BedTool(os.path.join(gtf_dir, 'proximal500_intron.bed'))
    }
    print("Genome annotation loaded.")

    # # process reference for use
    target = {
        "3UTR": ref_dict['3UTR'],
        "5UTR": ref_dict['5UTR'],
        "CDS": ref_dict['cds'],
        "other_exon": ref_dict['exon']-ref_dict['3UTR']-ref_dict['5UTR']-ref_dict['cds'],
        "px200_intron": ref_dict['proximal200'],
        "px500_intron": ref_dict['proximal500'].subtract(ref_dict['proximal200']),
        "distal_intron": ref_dict['intron'].subtract(ref_dict['exon']).subtract(ref_dict['proximal500'])
    }
    category_list = ['3UTR', '5UTR', 'CDS',
                     'other_exon', "px200_intron", "px500_intron", "distal_intron"]
    init = True

    print("Intersecting peaks with genome annotation...")
    for cat in category_list:
        bed_arr = []
        for interval in target[cat]:
            bed_arr.append('\t'.join([str(x) for x in interval.fields]))
            bed_arr[-1] = bed_arr[-1] + '\t' + cat
        bed_arr = list(dict.fromkeys(bed_arr))
        for i in range(len(bed_arr)):
            bed_arr[i] = bed_arr[i].split('\t')
        target[cat] = pybedtools.BedTool(bed_arr)

        if init:
            init = False
            result_bed = peaks.intersect(target[cat], wa=True, wb=True)
        else:
            result_bed = result_bed.cat(peaks.intersect(
                target[cat], wa=True, wb=True), postmerge=False)
    result_bed = result_bed.sort()

    print("Preparing output...")
    result_bed.saveas(outfn + '_')
    prepend = ['## Annotation peaks to genomic regions, all intersected genomic regions are presented.',
               '## CLAM version: %s'%config.__version__,
               '## Column 1:  Peak chromosome',
               '## Column 2:  Peak start',
               '## Column 3:  Peak end',
               '## Column 4:  Peak name',
               '## Column 5:  Peak score',
               '## Column 6:  Peak strand',
               '## Column 7:  Peak signal value',
               '## Column 8:  Peak pValue',
               '## Column 9:  Peak qValue',
               '## Column 10: Point-source called for this peak',
               '## Column 11: Genomic region chromosome',
               '## Column 12: Genomic region start',
               '## Column 13: Genomic region end',
               '## Column 14: Gene ID',
               '## Column 15: Quality score',
               '## Column 16: Genomic region strand',
               '## Column 17: Genomic region type']
    if os.path.exists(outfn):
        os.remove(outfn)
    for line in prepend:
        cmd = 'echo "{prepend}" >> {outfn}'.format(
            prepend=line, outfn=outfn)
        os.system(cmd)
    os.system('cat {outtmp} >> {outfn}'.format(
        outtmp=outfn + '_', outfn=outfn))
    os.remove(outfn+'_')
    print("DONE")


if __name__ == '__main__':
    # peak_fp, genome, outfn = sys.argv[1], sys.argv[2], sys.argv[3]
    os.chdir('/mnt/h/yi_lab/m6a/src/scripts/peakComposition')
    peak_in, genome, out_file = 'narrow_peak.unique.bed', 'mm10', 'annotate_peak.bed'
    if 'CLAM_DAT' not in os.environ or not download_data.check_genome_data(genome):
        print("Unable to find CLAM data folder for genomic regions, please try to download it using download_genome command.")
        print("Downloading...")
        download_data.download_genome(genome)
    genome_data = os.environ['CLAM_DAT']
    intersect_gtf_regions(
        peak_in, out_file, os.path.join(genome_data, genome))
