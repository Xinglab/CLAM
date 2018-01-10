#!/usr/bin/env python

"""
This re-aligner script is part of the CLAM pipeline.

It takes bam file as input, and outputs a weighed bam file for multi-mapped reads.

Tested under python 2.7.3
"""

__author__ = 'Zijun Zhang'
__version__ = '1.1.1'
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
from CLAM.preprocessor import *

logger = logging.getLogger('CLAM.Realigner')

class Bit:
	""" Binary Indexed Tree to store values in genomic intervals.
	Implementation modified from http://www.geeksforgeeks.org/binary-indexed-tree-or-fenwick-tree-2/
	Args:
	Returns:
	"""
	
	def __init__(self, n):
		sz = 1
		while n >= sz:
			sz *= 2
		self.size = sz
		self.array_size = n
		self.data = [0]*sz
		
	def sum(self, i):
		assert i >= 0
		if i==0:
			return 0
		if i > self.array_size:
			i = self.array_size 
		s = 0
		while i > 0:
			s += self.data[i]
			i -= i & -i
		return s
	
	def add(self, i, x):
		assert i > 0	
		while i < self.size:
			self.data[i] += x
			i += i & -i


def construct_BIT_track(subgraph, read_to_locations, ubam, unstranded=False):
	"""Construct BIT for each genomic region / node.
	Args:
		subgraph (list): a list of node names
		read_to_locations (dict): 
	Returns:
		node_track (Bit): a node-track dictionary
		multi_reads_weights (dict): a dictionary for multi-mapped reads.
	"""
	node_track = {}
	total_len = 0
	
	# initialized BIT tracks, add mreads to the tracks,
	# and keep a dict of read scores
	multi_reads_weights = defaultdict(dict)
	obs_reads = read_to_locations.keys()
	for read_x_qname in obs_reads:
		read_x_nodes = read_to_locations[read_x_qname]
		read_x_score = 1.0 / len(read_x_nodes)
		for node in read_x_nodes:
			chr, strand, start, end = node.split(':')
			start, end = int(start), int(end)
			if not node in node_track:
				this_len = end - start + 1
				node_track[node] = Bit(this_len)
			read_x_tag = read_x_nodes[node].opt('RT')
			node_locus = read_x_tag - start + 1
			node_track[node].add(node_locus, read_x_score)
			multi_reads_weights[read_x_qname][node]=[read_x_score, node_locus]
			#del read_to_locations[read_x_qname][node]
		#del read_to_locations[read_x_qname]
	
	# now add ureads by fetching from ubam; 
	# we don't need to keep track of them, just add the weights
	for node in node_track:
		chr, strand, start, end = node.split(':')
		start, end = int(start), int(end)
		is_reverse = True if strand=='-' else False
		uread_tagger = [x.opt('RT') for x in ubam.fetch(chr, start, end) \
			if unstranded or x.is_reverse==is_reverse]
		for uread_x_tagger in uread_tagger:
			if uread_x_tagger>=start and uread_x_tagger<=end:
				node_locus = uread_x_tagger - start + 1
				node_track[node].add(node_locus, 1)
	
	return node_track, multi_reads_weights



def run_EM(node_track, multi_reads_weights, w=50, epsilon=1e-6, max_iter=100, verbose=True):
	"""	EM implementation for re-assigning multi-mapped reads, given the 
	compatibility matrix of a subgraph.
	Args:
	Returns:
	"""
	iter=1
	residue=1
	while iter < max_iter and residue > epsilon:
		residue = 0
		reweight=defaultdict(dict)
		## calculate re-distribute probability; M-step
		for read in multi_reads_weights:
			for nd in multi_reads_weights[read]:
				track_len=node_track[nd].array_size
				old_score, read_tag = multi_reads_weights[read][nd]
				reweight[read][nd] = max( 0, node_track[nd].sum(min(track_len, read_tag + w)) - node_track[nd].sum(max(0,read_tag - w)) )
		## update track by expectation; E-step
		for read in reweight:
			dn=sum([reweight[read][x] for x in reweight[read]])
			if dn==0:
				logger.debug('Error: no read weight found @ %s.'%read )
				dn=1
			for nd in reweight[read]:
				old_score, read_tag = multi_reads_weights[read][nd]
				new_score = reweight[read][nd] / float(dn)
				node_track[nd].add(read_tag, new_score - old_score)
				residue += (old_score - new_score)**2
				multi_reads_weights[read][nd][0] = new_score
		if verbose and (not iter % 25 or iter == max_iter):
			logger.debug('Iter %d, residue = %f' % (iter, residue))
		iter += 1
	return multi_reads_weights


def build_read_cluster(alignment, chr_dict, location_to_reads, genomic_cluster_dict, unstranded=False, winsize=50):
	"""DOCSTRING
	Args:
	Returns:
	"""
	chr_list = chr_dict['name']
	chr_size = chr_dict['size']
	chrom = chr_list[alignment.reference_id]
	chr_len = chr_size[alignment.reference_id]
	site = alignment.opt('RT')
	is_reverse = alignment.is_reverse
	strand = '+' if unstranded or is_reverse==False else '-'
	this_mread_dict = {}
	this_mread_dict_set = defaultdict(set)
	discarded_mread_alignments = []
	
	## note to me: need to be more careful with
	## one read mapped to *multiple-locations* within one cluster
	## currently tossing away those alignments.. (in `discarded_mread_alignments`)
	
	## check junction reads; should be filtered out in tagging step
	if 'N' in alignment.cigarstring:
		return None, None, [alignment]
	
	## find the corresponding genomic cluster from `genomic_cluster_dict`
	chr_strand = chrom+':'+strand
	idx = bisect.bisect_right(genomic_cluster_dict[chr_strand], site)
	if not idx%2:
		print alignment
		raise Exception('%s falls out of region %s'%(alignment.qname, chr_strand+':'+str(site)) )
	start = genomic_cluster_dict[chr_strand][idx-1] - winsize
	start = 1 if start<1 else start
	end = genomic_cluster_dict[chr_strand][idx] + winsize
	end = chr_len-1 if end>=chr_len else end
	genomic_cluster = (chrom, strand, start, end)
	
	## fetch the reads
	# mread_list = [x for x in mbam.fetch(chrom, start, end) \
		# if (unstranded or x.is_reverse==is_reverse) and \
		# x.opt('RT')>=start and x.opt('RT')<=end]
	cluster_name = ':'.join([chrom, strand, str(genomic_cluster_dict[chr_strand][idx-1]), str(genomic_cluster_dict[chr_strand][idx])])
	mread_list = location_to_reads[cluster_name]
	for x in mread_list:
		this_mread_dict_set[x.qname].add(x)
	
	## find other mreads in this cluster
	for read_x_qname in this_mread_dict_set:
		if len(this_mread_dict_set[read_x_qname])>1:
			discarded_mread_alignments.extend( [ x for x in list(this_mread_dict_set[read_x_qname]) ])
		else:
			this_mread_dict[read_x_qname] = list(this_mread_dict_set[read_x_qname])[0]
	
	return genomic_cluster, this_mread_dict, discarded_mread_alignments


def construct_subgraph(location_to_reads, read_qname, mread_dict, processed_mreads, chr_dict, genomic_cluster_dict, winsize=50, unstranded=False):
	"""DOCSTRING
	Args:
		location_to_reads (dict): genomic cluster -> Alignment
		read_qname (str): target read ID
		mread_dict (dict): stores all read ID -> Alignment 
		processed_mreads (set): 
		chr_dict (dict): map ref_id to chrom_name, chrom_size
		genomic_cluster (dict): chrom:strand -> [interval1, interval2, ..]
	Returns:
	"""
	# record of processed alignments only need kept on within-subgraph level
	processed_mread_alignments = set()
	counter = 0
	# a list of `pysam.AlignedSegment` objects
	# note that all taggers are already stored in `pysam.AlignedSegment.opt('RT')`
	read_aln_list = [x for x in mread_dict[read_qname]] 
	processed_mreads.add(read_qname)
	read_to_locations = defaultdict(dict) # read_qname -> {node_name1:segment1, node_name2:segment2}
	
	# enumerate all connected components
	while True:
		counter+=1; #print "%i: %i"%(counter, len(read_aln_list))
		next_read_aln_list = []
		
		#gen = read_aln_list if len(read_aln_list)<200 else tqdm(read_aln_list)
		gen = read_aln_list
		for alignment in gen:
			## build a node for this mread alignment 
			## (if not already processed, i.e. built before)
			if alignment in processed_mread_alignments:
				continue
			
			genomic_cluster, this_mread_dict, discarded_mread_list = \
				build_read_cluster(alignment, chr_dict, 
					location_to_reads, genomic_cluster_dict, 
					unstranded=unstranded, winsize=winsize)
			_ = map(processed_mread_alignments.add, discarded_mread_list)
			if genomic_cluster is None:  # this cluster is invald (only double-mappers)
				continue
			
			## update location_to_reads, read2loc
			node_name = ':'.join([str(x) for x in genomic_cluster])
			#if node_name in subgraph:
				#logger.debug("I revisited '%s' at read '%s'."%(node_name, read_qname))
				#print("I revisited '%s' at read '%s'."%(node_name, read_qname))
				#break
			#subgraph.add(node_name)
			for x_qname in this_mread_dict:
				read_to_locations[x_qname].update({node_name :  this_mread_dict[x_qname]})
			
			## then add new alignments(edges) to generate connected nodes
			## in the next iteration
			_ = map(processed_mread_alignments.add, this_mread_dict.values())
			for read_x_qname in this_mread_dict:
				if read_x_qname in processed_mreads:
					continue
				x_aln_list = [aln for aln in mread_dict[read_x_qname] if not aln in processed_mread_alignments]
				next_read_aln_list.extend(x_aln_list)
			
			## .. and record to processed reads since we have generated
			## the nodes for them
			_ = map(processed_mreads.add, this_mread_dict.keys())
		
		# if no more connected nodes can be found, break loop 
		if len(next_read_aln_list)==0:
			break
		read_aln_list = next_read_aln_list		
	return read_to_locations, processed_mreads


def get_genomic_clusters(mbam, winsize=50, unstranded=False):
	"""Parsing the mbam to cluster the mread, and construct interval->alignment
	Args:
		mbam (pysam.Samfile): multi-read bam file handler
		winsize (int): 
	Returns:
	"""
	# chrom:+/- => [intv1_1, intv1_2, intv2_1, intv2_2]
	genomic_cluster_dict = defaultdict(list)
	# read_qname => [aln1, aln2, ..]
	mread_dict = defaultdict(list)
	# chrom:+/-:start:end => [read1_aln, read2_aln,]
	location_to_reads = defaultdict(list)
	chr_list=[x['SN'] for x in mbam.header['SQ']]
	chr_size=[x['LN'] for x in mbam.header['SQ']]
	chr_dict = {'name':chr_list, 'size':chr_size}
	logger.info('finding genomic clusters')
	for chrom in chr_list:
		## initialze the placeholder for current positive/negative strand clusters
		## pos/neg: [start, end, tag_counter]
		cur_cluster = {'+':[0,0,0], '-':[0,0,0]}
		cur_cluster_aln = {'+':[], '-':[]}
		for read_alignment in mbam.fetch(chrom):
			## should filter out junction reads in tagging step
			if 'N' in read_alignment.cigarstring:
				continue
			## add current alignment to mread_dict
			mread_dict[read_alignment.qname].append(read_alignment)
			site = read_alignment.opt('RT')
			strand = '-' if read_alignment.is_reverse and not unstranded else '+'
			### if this read is within the window size
			if site <= cur_cluster[strand][1]+winsize:
				if site < cur_cluster[strand][0]:
					cur_cluster[strand][0] = site
				if site > cur_cluster[strand][1]:
					cur_cluster[strand][1] = site
				cur_cluster[strand][2] += 1
				cur_cluster_aln[strand].append(read_alignment)
			### otherwise, push the current cluster to `genomic_cluster_dict`
			else:
				if cur_cluster[strand][2] > 0:
					genomic_cluster_dict[chrom+':'+strand].extend(
						[
							cur_cluster[strand][0], 
							cur_cluster[strand][1]+1
						]
						)
					genomic_cluster_str = ':'.join([chrom, strand, str(cur_cluster[strand][0]), str(cur_cluster[strand][1]+1) ])
					location_to_reads[genomic_cluster_str].extend(cur_cluster_aln[strand])
				cur_cluster[strand] = [site, site+1, 1]
				cur_cluster_aln[strand] = [read_alignment]
		## remember to push the last genomic cluster to dict
		if cur_cluster['+'][2] > 0:
			genomic_cluster_dict[chrom+':+'].extend([cur_cluster['+'][0],cur_cluster['+'][1]+1])
			genomic_cluster_str = ':'.join([chrom, '+', str(cur_cluster['+'][0]), str(cur_cluster['+'][1]+1) ])
			location_to_reads[genomic_cluster_str].extend(cur_cluster_aln['+'])
		if cur_cluster['-'][2] > 0:
			genomic_cluster_dict[chrom+':-'].extend([cur_cluster['-'][0],cur_cluster['-'][1]+1])
			genomic_cluster_str = ':'.join([chrom, '-', str(cur_cluster['-'][0]), str(cur_cluster['-'][1]+1) ])
			location_to_reads[genomic_cluster_str].extend(cur_cluster_aln['-'])
	
	return genomic_cluster_dict, mread_dict, location_to_reads



def realigner(in_bam, out_dir, max_hits=100, max_tags=-1, read_tagger_method='median', 
		winsize=50, unstranded=False, retag=False):
	"""DOCSTRING
	Args:
	Returns:
	"""
	# logging the parameter values
	frame = inspect.currentframe()
	args, _, _, values = inspect.getargvalues(frame)
	msg = 'Params:\n'
	for i in args:
		msg += "%s = %s \n"%(i, values[i])
	logger.info(msg)
	# preprocessing
	if retag or not (
			os.path.isfile(os.path.join(out_dir,'unique.sorted.bam')) and \
			os.path.isfile(os.path.join(out_dir,'multi.sorted.bam')) \
			) :
		filter_bam_multihits(in_bam, max_tags=max_tags, max_hits=max_hits, out_dir=out_dir, read_tagger_method=read_tagger_method, 
			omit_detail=True)

	# file handlers
	if max_tags>0:
		mbam = pysam.Samfile(os.path.join(out_dir, 'multi.sorted.collapsed.bam'),'rb')
		ubam = pysam.Samfile(os.path.join(out_dir, 'unique.sorted.collapsed.bam'),'rb')
	else:
		mbam = pysam.Samfile(os.path.join(out_dir, 'multi.sorted.bam'),'rb')
		ubam = pysam.Samfile(os.path.join(out_dir, 'unique.sorted.bam'),'rb')
	obam = pysam.Samfile(os.path.join(out_dir, 'realigned.bam'), 'wb', template = mbam)
	chr_list=[x['SN'] for x in ubam.header['SQ']]
	chr_size=[x['LN'] for x in mbam.header['SQ']]
	chr_dict = {'name':chr_list, 'size':chr_size}
	
	# construct the `mread_dict`, this will be needed throughout the program;
	# also construct the genomic cluster dict and cluster to alignment,
	# by going through all mreads at once
	genomic_cluster_dict, mread_dict, location_to_reads = get_genomic_clusters(mbam, winsize=winsize, unstranded=unstranded)
	
	# keep a record of processed reads
	processed_mreads = set()
	
	# iterate through all mreads
	logger.info('running em')
	for read_qname in mread_dict:
		if read_qname in processed_mreads:
			continue
			
		## construct the fully-connected subgraph for each read
		read_to_locations, processed_mreads = \
			construct_subgraph(location_to_reads, read_qname, mread_dict, processed_mreads, chr_dict, \
				genomic_cluster_dict, winsize=winsize, unstranded=unstranded)
		subgraph = set()
		for read in read_to_locations:
			_ = map(subgraph.add, read_to_locations[read].keys())
		subgraph = list(subgraph)
		logger.debug("|v|=%i, |e|=%i"%(len(subgraph), len(read_to_locations)) )
		
		## build the BIT tracks
		node_track, multi_reads_weights = \
			construct_BIT_track(subgraph, read_to_locations, ubam, unstranded)
			
		## run EM
		multi_reads_weights = \
			run_EM(node_track, multi_reads_weights, w=winsize)
		
		## write to obam
		for read in multi_reads_weights:
			for node in multi_reads_weights[read]:
				alignment = read_to_locations[read][node]
				score = round(multi_reads_weights[read][node][0], 3)
				alignment.set_tag('AS', score)
				#alignment.set_tag('PG', 'CLAM')
				obam.write(alignment)
	# sort the final output
	logger.info('sorting output')
	obam.close()
	ubam.close()
	mbam.close()
	obam_sorted_fn = os.path.join(out_dir, 'realigned.sorted.bam')
	pysam.sort('-o', obam_sorted_fn, os.path.join(out_dir, 'realigned.bam'))
	pysam.index(obam_sorted_fn)
	os.remove(os.path.join(out_dir, 'realigned.bam'))
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
		max_tags = args.max_tags
		retag = args.retag
		winsize = args.winsize
		unstranded = args.unstranded
		
		#logger = logging.getLogger('CLAM.Realigner')
		logger.info('start')
		logger.info('run info: %s'%(' '.join(sys.argv)))
		
		realigner(in_bam, out_dir, max_hits=max_hits, max_tags=max_tags, read_tagger_method=tag_method, 
			winsize=winsize, unstranded=unstranded, retag=retag)
		
		logger.info('end')
	except KeyboardInterrupt():
		sys.exit(0)
	return
	
	
	

if __name__=='__main__':		
	### set up logger
	logger = logging.getLogger('CLAM')
	logger.setLevel(logging.DEBUG)
	# create file handler which logs even debug messages
	fh = logging.FileHandler(
		'CLAM.Realigner.'+'-'.join(str(datetime.datetime.now()).replace(':','-').split()) + '.log')
	fh.setLevel(logging.INFO)
	# create console handler with a higher log level
	ch = logging.StreamHandler()
	ch.setLevel(logging.DEBUG)
	# create formatter and add it to the handlers
	formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s -\n %(message)s')
	fh.setFormatter(formatter)
	ch.setFormatter(formatter)
	# add the handlers to the logger
	logger.addHandler(fh)
	logger.addHandler(ch)
	logger.info('start')
	
	logger.info('run info: %s'%(' '.join(sys.argv)))
	bam, out_dir, out_dir = sys.argv[1:4]
	retag = False
	if len(sys.argv)>4:
		tagger_method = sys.argv[4]
		retag = True
		logger.info('retagging with "%s"'%tagger_method)
	else:
		tagger_method = 'median'
		
	if retag or not (
			os.path.isfile(os.path.join(out_dir,'unique.sorted.bam')) and \
			os.path.isfile(os.path.join(out_dir,'multi.sorted.bam')) \
			) :
		filter_bam_multihits(bam, max_hits=100, out_dir=out_dir, read_tagger=lambda x: read_tagger(x, tagger_method))
	realigner(out_dir, out_dir, winsize=50, unstranded=False)
	logger.info('end')