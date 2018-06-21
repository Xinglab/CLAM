#!/usr/bin/env python

"""This re-aligner script is part of the CLAM pipeline.

This subcommand will run expectation-maxmization to assign the multi-mapped reads in a probablistic framework. 
More details about the EM model is described in our NAR paper.

Note when `--retag` is specified, `realigner` will re-run `preprocessor` regardless; otherwise, it will use 
the prepared files in `outdir` if available.

Example run:
	```
	CLAM realigner -i path/to/input/Aligned.out.bam -o path/to/clam/outdir/ --read-tagger-method start --retag
	```

Tested under python 2.7
"""

__author__ = 'Zijun Zhang'
__version__ = '1.1.2'
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
from .preprocessor import *

logger = logging.getLogger('CLAM.Realigner')

class Bit:
	""" Binary Indexed Tree to store values in genomic intervals.
	Implementation modified from http://www.geeksforgeeks.org/binary-indexed-tree-or-fenwick-tree-2/
	Args:
		n (int): length of the interval to construct
	Returns:
		a BIT object with `add` and `sum` method over arbitrary sub-intervals
		with O(log(n)) time
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
		node_track (dict/BIT): a node-track dictionary; node_name => BIT
		multi_reads_weights (dict): a dictionary for multi-mapped reads; read_qname => node => [score, locus]
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
		node_track (dict): dict. of BIT returned from `construct_BIT_track`
		multi_reads_weights (dict): dict of mread qname and locus returned from `construct_BIT_track`
		w (int): window size for search vicinity reads
		epsilon (float): a small number for testing convergence between iterations
		max_iter (int): maximum iterations of EM
		verbose (bool, options): prints status in verbose mode
	Returns:
		multi_reads_weights (dict): the mread weight after EM
	"""
	iter = 1
	residue = 1
	#n_est = sum([1. for r in multi_reads_weights for n in multi_reads_weights[r] ])
	n_est = sum([1. for r in multi_reads_weights ])
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
			residue /= n_est
		if verbose and (not iter % 10 or iter == max_iter):
			logger.debug('Iter %d, residue = %f' % (iter, residue))
		iter += 1
	return multi_reads_weights


def build_read_cluster(alignment, chr_dict, location_to_reads, genomic_cluster_dict, unstranded=False, winsize=50):
	"""Given an alignment, find its genomic cluster, and all other mreads 
	in that cluster
	Args:
		alignment (pysam.AlignedSegment): pysam alignment object
		chr_dict (dict): a dict of chrom name and sizes
		location_to_reads (dict): stores all mreads indexed by aligned locus; cluster name => mread alignments
		genomic_cluster_dict (dict): stores genomic clusters; chrom => [intv1, intv2, ..]
		unstranded (bool): if true, don't use the strand info in alignment
		winsize (int): window size for search ureads
	Returns:
		genomic_cluster (tuple): the target chrom and coordinates after expanding the window size
		this_mread_dict (dict): dict. of mread alignments in the target cluster indexed by read_qname
		discarded_mread_alignments (list): discarded mread alignments because of multiple occurences within one cluster
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
	#if 'N' in alignment.cigarstring:
	#	return None, None, [alignment]
	
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
	cluster_name = ':'.join([chrom, strand, str(genomic_cluster_dict[chr_strand][idx-1]), str(genomic_cluster_dict[chr_strand][idx])])
	if not cluster_name in location_to_reads:
		raise Exception("cannot find cluster %s in `location_to_reads`"%cluster_name)
	mread_list = location_to_reads[cluster_name]
	del location_to_reads[cluster_name]
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
	"""Given a mread_qname, find exhaustively all other connected mreads.
	Args:
		location_to_reads (dict): genomic cluster => Alignment
		read_qname (str): target read ID
		mread_dict (dict): stores all read ID => Alignment 
		processed_mreads (set): 
		chr_dict (dict): map ref_id to chrom_name, chrom_size
		genomic_cluster (dict): chrom:strand => [interval1, interval2, ..]
	Returns:
		read_to_locations (dict): collect a subset of mread alignments in the same 
			subgraph starting with read_qname 
		processed_mreads (set): record all processed mread_qname to avoid re-processing
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
			## record those discarded alignments/reads
			## note: we mark discarded_mread as processed as well,
			## so as not to create a bias to less clustered regions.
			_ = map(processed_mread_alignments.add, discarded_mread_list)
			_ = map(processed_mreads.add, [x.qname for x in discarded_mread_list])
			if genomic_cluster is None:  # this cluster is invald (only double-mappers)
				continue
			
			## update read_to_locations
			node_name = ':'.join([str(x) for x in genomic_cluster])
			#if node_name in subgraph:
				#logger.debug("I revisited '%s' at read '%s'."%(node_name, read_qname))
				#print("I revisited '%s' at read '%s'."%(node_name, read_qname))
				#break
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
	"""Parsing the mbam to cluster the mread, and construct interval=>alignment.
	Using the same object in difference references, and just keep one copy of 
	the mread-alignments to minimize memory usage.
	Args:
		mbam (pysam.Samfile): multi-read bam file handler
		winsize (int): window size for search mreads
		unstranded (bool): if turned on, all reads will be pushed to forward strand
	Returns:
		genomic_cluster_dict (dict): chrom:+/- => [intv1, intv2, ..]
		mread_dict (dict): read_qname => [aln1, aln2, ..]
		location_to_reads (dict): chrom:strand:start:end => [read1_aln, real2_aln, ..]
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
			#if 'N' in read_alignment.cigarstring:
			#	continue
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
	"""The main entry for CLAM-realigner.
	Args:
		in_bam (str): filepath for input bam
		out_dir (str): filepath for CLAM output folder
		max_hits (int):  maximum number of aligned loci allowed for mreads
		max_tags (int): maximum number of identical alignments allowed for each 
			genomic locus, more amount will be collapsed; -1 is no collapsing
		read_tagger_method (str): the tagger function type
		winsize (int): window size 
		unstranded (bool): ignore alignment strand info if turned on
		retag (bool): force to call `preprocessor` to process `in_bam` if turned on
	Returns:
		None
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
	else:
		logger.info("found existing bams; skipped tagging.")

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
	logger.debug('found %i mreads @ %i locations' % ( len(mread_dict), len(location_to_reads) ) )
	
	# keep a record of processed reads
	processed_mreads = set()
	
	# iterate through all mreads
	logger.info('running em')
	subg_counter = 0
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
		#if len(subgraph)==1 and len(read_to_locations)>10:
		#	raise Exception('Incorrect mread assigned to one location')
		if len(subgraph)==0:
			continue
		subg_counter += 1
		logger.debug("subgraph %i: |e|=%i, |v|=%i"%(subg_counter, len(read_to_locations), len(subgraph)) )
		
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
	"""The command-line parser for CLAM-realigner
	Args:
		args (argparse.ArgumentParser): receives commandline arguments
	Returns:
		None
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
	# *****
	# NOTE: below is used for debugging purpose;
	# users should call from `CLAM subcommand` instead
	# of running this script directly
	# *****
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