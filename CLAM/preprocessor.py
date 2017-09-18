#!/usr/bin/env python

"""
This preprocessing script is part of the CLAM pipeline.

It takes bam file as input, separates unique- and multi-mapped reads, and
tag the read to a specific locus given a tagging method.

Tested under python 2.7.3
"""

__author__ = 'Zijun Zhang'
__version__ = '1.1.0'
__email__ = 'zj.z@ucla.edu'

import os
import sys
import pysam
import numpy as np
from collections import defaultdict
from tqdm import tqdm
import logging
import datetime
import bisect
import argparse as ap
import inspect

logger = logging.getLogger('CLAM.Preprocessor')

def read_tagger(alignment, method='median'):
	""" tag a read alignment to a genomic locus
	Args:
	Returns:
	"""
	tagger_func = {
		# center of the read; must dicard junction reads
		'median': lambda x: -1 if 'N' in x.cigarstring else int(np.median(x.positions))+1,
		# start site of the read; trunction in iCLIP/eCLIP
		'start': lambda x: x.positions[-1] if x.is_reverse else x.positions[0]+1
		}
	try:
		tag=tagger_func[method](alignment)
	except:
		tag=-1
	return tag



def filter_bam_multihits(filename, max_hits, out_dir, read_tagger, omit_detail=True):
	"""Pre-processing function for cleaning up the input bam file.
	Args:
	Returns:
	"""
	logger.info('filtering input bam')
	
	in_bam = pysam.Samfile(filename,'rb')
	# unique read bam
	ubam_fn = os.path.join(out_dir, 'unique.bam')
	sorted_ubam_fn = os.path.join(out_dir, 'unique.sorted.bam')
	ubam=pysam.Samfile(ubam_fn, 'wb', template=in_bam)
	unique_counter = 0
	
	# multi-read bam
	mbam_fn = os.path.join(out_dir, 'multi.bam')
	sorted_mbam_fn = os.path.join(out_dir, 'multi.sorted.bam')
	mbam=pysam.Samfile(mbam_fn, 'wb', template=in_bam)
	mread_set = set()
	
	# splitting unique and multi- reads
	# and add the read taggers we need
	#for read in tqdm(in_bam):
	counter = 0
	for read in in_bam:
		# poor man's progress bar
		counter += 1
		if not counter % 10**6:
			logger.debug('tagged %i alignments'%counter)
		read_tag = read_tagger(read)
		## skip reads with unassigned tagger
		if read_tag==-1:
			continue
		read.tags += [('RT', read_tag)] ## add the tag
		
		## omit the details in read sequence and quality
		## recommended for larger bam because this
		## can save some memory/storage for large bams
		if omit_detail:
			read.query_sequence = '*'
			read.query_qualities = '0'
		
		if read.is_secondary or (read.has_tag('NH') and read.opt("NH")>1):
			try:
				if read.opt("NH") < max_hits:
					mbam.write(read)
					mread_set.add(read.qname)
			except KeyError:
				#print read
				raise Exception('%s: missing NH tag when is_secondary=%s'%(read.qname,read.is_secondary))
		else:
			ubam.write(read)
			unique_counter += 1
	
	in_bam.close()
	ubam.close()
	mbam.close()
	
	# sorting
	pysam.sort('-o', sorted_ubam_fn, ubam_fn)
	os.remove(ubam_fn)
	pysam.sort('-o', sorted_mbam_fn, mbam_fn)
	os.remove(mbam_fn)
	pysam.index(sorted_ubam_fn)
	pysam.index(sorted_mbam_fn)
	
	# log the statistics
	multi_counter = len(mread_set)
	logger.info(
			'Unique reads = %s;  ' % unique_counter + \
			'Multi reads = %s (%.2f %%)' % \
			( multi_counter, float(multi_counter)/(multi_counter+unique_counter)*100 )
		)
	return


def parser(args):
	"""DOCSTRING
	Args
	Returns
	"""
	try:
		in_bam = args.in_bam
		out_dir  = args.out_dir
		if not os.path.isdir(out_dir):
			os.mkdir(out_dir)
		tag_method = args.tag_method
		max_hits = args.max_hits
		
		#logger = logging.getLogger('CLAM.Preprocessor')
		logger.info('start')
		logger.info('run info: %s'%(' '.join(sys.argv)))
		
		filter_bam_multihits(in_bam, max_hits=max_hits, out_dir=out_dir, read_tagger=lambda x: read_tagger(x, method=tag_method))
		
		logger.info('end')
	except KeyboardInterrupt():
		sys.exit(0)
	return