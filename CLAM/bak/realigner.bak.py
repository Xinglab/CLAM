#!/usr/bin/env python

"""
This re-aligner script is part of the CLAM pipeline.

It takes bam file as input, and outputs a weighed bam file for multi-mapped reads.

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

### get logger
###
logger = logging.getLogger('CLAM.Realigner')
###


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
	Returns:
	Returns a node-track dictionary and a dictionary for multi-mapped reads.
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
	
	# now add ureads by fetching from ubam
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
		if verbose and (not iter % 10 or iter == max_iter):
			logger.debug('Iter %d, residue = %f' % (iter, residue))
		iter += 1
	return multi_reads_weights


def build_read_cluster(alignment, chr_list, mbam, unstranded=False, winsize=50):
	"""DOCSTRING
	Args:
	Returns:
	"""
	chrom = chr_list[alignment.reference_id]
	site = alignment.opt('RT')
	is_reverse = alignment.is_reverse
	this_mread_dict = {}
	this_mread_dict_set = defaultdict(set)
	discarded_mread_alignments = []
	## note to me: need to be more careful with
	## one read mapped to *multiple-locations* within one cluster
	## currently tossing away those alignments.. (in `discarded_mread_alignments`)
	
	# find the right boundary
	current = site
	while True:
		mread_list = [x for x in mbam.fetch(chrom, current, current+winsize) \
			if unstranded or x.is_reverse==is_reverse]
		for x in mread_list:
			this_mread_dict_set[x.qname].add(x)
		mread_tagger = [x.opt('RT') for x in mread_list]
		if len(mread_tagger)==0:
			break
			#return (None,None,discarded_mread_alignments)
		end = max(mread_tagger) + winsize
		if max(mread_tagger)!=current:
			## has to restrict step-size to smaller than winsize; 
			## in order to avoid missing of very short reads
			current = max(mread_tagger) if max(mread_tagger)<current+winsize else current+winsize
		else:
			break
		
	# find the left boundary
	current = site
	while True:
		mread_list = [x for x in mbam.fetch(chrom, current-winsize, current) \
			if unstranded or x.is_reverse==is_reverse]
		for x in mread_list:
			this_mread_dict_set[x.qname].add(x)
		mread_tagger = [x.opt('RT') for x in mread_list]
		if len(mread_list)==0:
			break
			#return (None,None,discarded_mread_alignments)
		start = min(mread_tagger) - winsize
		if min(mread_tagger) != current:
			current = min(mread_tagger) if min(mread_tagger)>current-winsize else current-winsize
		else:
			break
		
	strand = '+' if unstranded or is_reverse==False else '-'
	for read_x_qname in this_mread_dict_set:
		if len(this_mread_dict_set[read_x_qname])>1:
			discarded_mread_alignments.extend( [ x for x in list(this_mread_dict_set[read_x_qname]) ])
		else:
			this_mread_dict[read_x_qname] = list(this_mread_dict_set[read_x_qname])[0]
		
	genomic_cluster = (chrom, strand, start, end)
	
	return genomic_cluster, this_mread_dict, discarded_mread_alignments


def construct_subgraph(mbam, read_qname, mread_dict, processed_mreads, chr_list, winsize=50, unstranded=False):
	"""DOCSTRING
	Args:
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
		counter+=1; print "%i: %i"%(counter, len(read_aln_list))
		next_read_aln_list = []
		
		gen = read_aln_list if len(read_aln_list)<200 else tqdm(read_aln_list)
		for alignment in gen:
			## build a node for this mread alignment 
			## (if not already processed, i.e. built before)
			if alignment in processed_mread_alignments:
				continue
			
			genomic_cluster, this_mread_dict, discarded_mread_list = \
				build_read_cluster(alignment, chr_list, mbam, unstranded=unstranded, winsize=winsize)
			_ = map(processed_mread_alignments.add, discarded_mread_list)
			if genomic_cluster is None:  # this cluster is invald (only double-mappers)
				continue
			
			## update loc2read, read2loc
			node_name = ':'.join([str(x) for x in genomic_cluster])
			#if node_name in subgraph:
			#	logger.debug("I revisited '%s' at read '%s'."%(node_name, read_qname))
			#	break
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


def realigner(out_dir, tmp_dir, winsize=50, unstranded=False):
	"""DOCSTRING
	Args:
	Returns:
	"""
	# file handlers
	mbam = pysam.Samfile(os.path.join(tmp_dir, 'multi.sorted.bam'),'rb')
	ubam = pysam.Samfile(os.path.join(tmp_dir, 'unique.sorted.bam'),'rb')
	obam = pysam.Samfile(os.path.join(out_dir, 'realigned.bam'), 'wb', template = mbam)
	chr_list=[x['SN'] for x in ubam.header['SQ']]
	
	# construct the mread_dict; this will be needed throughout
	mread_dict = defaultdict(list)
	for alignment in mbam:
		mread_dict[alignment.qname].append(alignment)
		
	# keep a record of processed reads
	processed_mreads = set()
	
	# iterate through all mreads
	for read_qname in mread_dict:
		if read_qname in processed_mreads:
			continue
			
		## construct the fully-connected subgraph for each read
		read_to_locations, processed_mreads = \
			construct_subgraph(mbam, read_qname, mread_dict, processed_mreads, chr_list, winsize=winsize, unstranded=unstranded)
		subgraph = set()
		for read in read_to_locations:
			_ = map(subgraph.add, read_to_locations[read].keys())
		subgraph = list(subgraph)
		
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
				alignment.set_tag('PG', 'CLAM')
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
	
def read_tagger(alignment, method='median'):
	""" tag a read alignment to a genomic locus
	Args:
	Returns:
	"""
	tagger_func = {
		'median': lambda x:  int(np.median(x.positions))+1,
		'start': lambda x: x.positions[-1] if x.is_reverse else x.positions[0]+1
		}
	try:
		tag=tagger_func[method](alignment)
	except:
		tag=-1
	return tag


def filter_bam_multihits(filename, max_hits, tmp_dir, read_tagger, omit_detail=True):
	"""Pre-processing function for cleaning up the input bam file.
	Args:
	Returns:
	"""
	logger.info('Filtering input bam..')
	
	in_bam = pysam.Samfile(filename,'rb')
	# unique read bam
	ubam_fn = os.path.join(tmp_dir, 'unique.bam')
	sorted_ubam_fn = os.path.join(tmp_dir, 'unique.sorted.bam')
	ubam=pysam.Samfile(ubam_fn, 'wb', template=in_bam)
	unique_counter = 0
	
	# multi-read bam
	mbam_fn = os.path.join(tmp_dir, 'multi.bam')
	sorted_mbam_fn = os.path.join(tmp_dir, 'multi.sorted.bam')
	mbam=pysam.Samfile(mbam_fn, 'wb', template=in_bam)
	mread_set = set()
	
	# splitting unique and multi- reads
	# and add the read taggers we need
	for read in tqdm(in_bam):
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
			read.query_qualities = [0]
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
	bam, tmp_dir, out_dir = sys.argv[1:4]
	retag = False
	if len(sys.argv)>4:
		tagger_method = sys.argv[4]
		retag = True
		logger.info('Retag with "%s"'%tagger_method)
	else:
		tagger_method = 'median'
		
	if retag or not (
			os.path.isfile(os.path.join(tmp_dir,'unique.sorted.bam')) and \
			os.path.isfile(os.path.join(tmp_dir,'multi.sorted.bam')) \
			) :
		filter_bam_multihits(bam, max_hits=100, tmp_dir=tmp_dir, read_tagger=lambda x: read_tagger(x, tagger_method))
	realigner(out_dir, tmp_dir, winsize=50, unstranded=False)
	logger.info('end')