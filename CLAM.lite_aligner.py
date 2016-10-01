#!/usr/bin/python

"""
This re-aligner script is part of the CLAM pipeline.

It takes bam file as input, and outputs a weighed bam file for multi-mapped reads.

Tested under python 2.7.3
"""

__author__ = 'Zijun Zhang'
__version__ = '1.0.0'
__email__ = 'zj.z@ucla.edu'


from optparse import OptionParser
import os, subprocess, shutil, sys, copy
import pysam
import numpy as np
from collections import defaultdict, deque
from time import strftime
import bisect, operator
import cPickle as pickle
import pybedtools


def main():
	"""
	The main wrapper for CLAM re-aligner.
	"""
	# options parsing
	usage='usage: %prog <options> input_file.bam'
	parser=OptionParser(usage)
	
	parser.add_option('-o', dest='output_dir', default='./out_CLAM', help='Output file folder [Default %default]')
	parser.add_option('-t', dest='tmp_dir', default='./tmp_CLAM', help='Temporary file folder [Default %default]')
	parser.add_option('-w', dest='window_size', type='int', default=50, help='Local window size for EM [Default: %default]')
	parser.add_option('--max-multihits', dest='max_multihits', type='int', default=100, help='Discard reads mapped to more than <max_multihits> locations. [Default: %default]')
	parser.add_option('--min-unique-reads', dest='min_unique_reads', type='int', default=0, help='Discard genomic regions with less than <min_unique_reads> of unique reads. [Default: %default]')
	parser.add_option('--is-stranded', dest='is_stranded', default=False, action='store_true', help='Indicates if the reads are mapped with strand information. [Default: %default]')
	parser.add_option('--resume', dest='resume', action='store_true', default=False, help='Resume mode - skipping pre-processing [Default: %default]')
	parser.add_option('--verbose', dest='verbose', action='store_true', default=False, help='Verbose mode - print out all intermediate steps [Default: %default]')
	parser.add_option('--max-gap', dest='max_gaps', type='int', default=50, help='Maximum distance allowed in grouping reads. [Default: %default]')
	#parser.add_option('-g', dest='gtf', default='/u/home/f/frankwoe/scratch/hg19_gencodeV19.sorted.bed', help='GTF file (only need to specified when using -c). [Default: %default]')
	#parser.add_option('--covariate-site-min', dest='cov_site_min', type='float', default=0, help='Minumim value required for same genomic region in covariate file [Default: %default]')
	#parser.add_option('--covariate-gene-min', dest='cov_gene_min', type='float', default=1.0, help='Minumim value required for same gene in covariates (e.g. RPKM/FPKM). [Default: %default]')
	
	(options,args)=parser.parse_args()
	
	if len(args)<1:
		parser.error('Missing required input.')
	
	input_file=args[0]
	input_cov = args[1] if len(args)>1 else None
	
	output_dir=os.path.abspath(options.output_dir)
	tmp_dir=os.path.abspath(options.tmp_dir)
	verbose=options.verbose
	steps_finished=[]
	
	print_time_stamp('Program start.')
	# checking for output files that already exist; if non-exist, call preprocessing subroutines
	# if resume is turned on, start from finished files
	# NOTE: currently doesn't check if a file is created but truncated
	if options.resume:
		if os.path.isfile('%s/filter%d.sorted.bam' % (tmp_dir, options.max_multihits)):
			filter_filename='%s/filter%d.sorted.bam' % (tmp_dir,options.max_multihits)
			multiread_set=pickle.load(open(tmp_dir + '/multiread_set.pdata','rb'))
			steps_finished.append('filtered%d.sorted.bam' % options.max_multihits)
		if os.path.isfile(tmp_dir + '/genomic_regions_%s_%s.pdata' % (str(options.max_gaps), str(options.min_unique_reads))):
			genomic_regions=pickle.load(open(tmp_dir + '/genomic_regions_%s_%s.pdata' % (str(options.max_gaps), str(options.min_unique_reads) ),'rb'))
			location_to_reads=pickle.load(open(tmp_dir + '/loc2read.pdata','rb'))
			read_to_locations=pickle.load(open(tmp_dir + '/read2loc.pdata','rb'))
			steps_finished.extend(['genomic_regions.pdata', 'loc2read.pdata', 'read2loc.pdata'])
		#if os.path.isfile(tmp_dir + '/preserved_nodes_%s_%s.pdata' % (str(options.cov_site_min), str(options.cov_gene_min)) ) and not input_cov is None:
		#	preserved_nodes=pickle.load(open(tmp_dir + '/preserved_nodes_%s_%s.pdata' % (str(options.cov_site_min), str(options.cov_gene_min)),'rb'))
		#	steps_finished.append('preserved_nodes.pdata')
		if os.path.isfile(output_dir + '/CLAM_mapper.out'):
			steps_finished.append('CLAM_mapper.out')
		print_time_stamp('Resume mode On, found files: ' + ','.join(steps_finished))
	else:
		if(os.path.isdir(tmp_dir)):
			shutil.rmtree(tmp_dir)
		if(os.path.isdir(output_dir)):
			shutil.rmtree(output_dir)
		os.mkdir(tmp_dir)
		os.mkdir(output_dir)
		
	write_parameter_log(options, args, output_dir)
	if not ( options.resume and 'filter_filename' in locals() ):
		filter_filename, multiread_set=filter_multihits(input_file, options.max_multihits, verbose, tmp_dir)
	
	bamfile=pysam.Samfile(filter_filename,'rb')
	
	# Construct genomic regions and read-location / location-read dictionary
	if not ( options.resume and 'genomic_regions.pdata' in steps_finished ):
		genomic_regions, location_to_reads, read_to_locations=get_genomic_regions(bamfile, options.max_gaps, verbose, tmp_dir, options.is_stranded, options.min_unique_reads, multiread_set)
	
	
	bamfile.close()	
	if not input_cov is None:
		#if  not (options.resume and 'preserved_nodes.pdata' in steps_finished ):
		#	preserved_nodes = filter_by_covr(input_cov, options.cov_site_min, options.cov_gene_min, options.gtf, genomic_regions, tmp_dir)
		#	pickle.dump(preserved_nodes, open(tmp_dir+'/preserved_nodes_%s_%s.pdata' % (str(options.cov_site_min), str(options.cov_gene_min)), 'wb'), -1)
		genomic_regions, read_to_locations, location_to_reads = update_by_filter(genomic_regions, read_to_locations, location_to_reads, preserved_nodes)
	
	print_time_stamp('Pre-process done.')
	
	if options.resume and 'CLAM_mapper.out' in steps_finished:
		nodes_finished = read_rm_out(output_dir+'/CLAM_mapper.out')
	else:
		nodes_finished=set()
	
	# Call EM model to assign multi-mapped reads
	print_time_stamp('EM start.')
	iter=0
	out=open(output_dir + '/CLAM_mapper.out','w',0) if options.resume else open(output_dir + '/CLAM_mapper.out','a',0)
	seen=nodes_finished
	for chr_strand in genomic_regions:
		chr, strand = chr_strand.split(':')
		for start, end in genomic_regions[chr_strand]:
			node=chr_strand + ':' + str(start) + ':' + str(end)
			k = len(seen)
			if not (k+1) % 1000:
					print_time_stamp(str(k+1) + ' finished.')
			subgraph, seen=search_node_subg(node, location_to_reads, read_to_locations, seen)
			if subgraph is None or len(subgraph)<2:
				continue
			
			node_track, multi_reads_weights = construct_track_lite(subgraph, location_to_reads, read_to_locations)

			iter += 1
			if len(multi_reads_weights)<1:
				print_time_stamp('Error occured for: ' + ','.join(subgraph) + ': No fully contained reads found.' + str(len(obs_reads)))
				continue
			if verbose:
				print_time_stamp('Round ' + str(iter) + ': seen = ' + str(len(seen)) + '; current subgraph = ' + str(len(subgraph)) + '; obs reads = ' + str(len(multi_reads_weights)))

			
			new_reads_weights = runEM(node_track, multi_reads_weights, w=options.window_size)
			wrt_content = make_write_content(new_reads_weights)
			out.write(wrt_content)
			
	out.close()
	# write output files
	print_time_stamp('Sorting output Bedfile.')
	subprocess.call(''' sort -k1,1 -k2,2n %s/CLAM_mapper.out > %s/CLAM_mapper.sorted.out ''' % (output_dir, output_dir), shell=True)
	header_cmd='samtools view -H ' + tmp_dir + '/filter100.sorted.bam > ' + output_dir + '/sam_header.sam'
	subprocess.call(header_cmd, shell=True)
	body_cmd = ''' awk '{if($6=="+"){print $4"\t256\t"$1"\t"$2+1"\t0\t"$3-$2+1"M\t*\t0\t0\t*\t*\tAS:f:"$5}else{print $4"\t272\t"$1"\t"$2+1"\t0\t"$3-$2+1"M\t*\t0\t0\t*\t*\tAS:f:"$5 }}' ''' + output_dir + '/CLAM_mapper.sorted.out > ' + output_dir + '/CLAM_mapper.sorted.sam'
	subprocess.call(body_cmd, shell=True)
	makeBam_cmd = 'cat %s/sam_header.sam %s/CLAM_mapper.sorted.sam | samtools view -bS - > %s/assigned_multimapped_reads.bam' % (output_dir, output_dir,output_dir)
	subprocess.call(makeBam_cmd, shell=True)
	index_cmd = 'samtools index %s/assigned_multimapped_reads.bam' % output_dir
	subprocess.call(index_cmd, shell=True)
	print_time_stamp('Re-alignment is done.')

def write_parameter_log(options, args, output_dir):
	"""
	Write paramter values to a log file, named by current time.
	"""
	with open(output_dir+'/CLAM_Aligner.Log.'+ strftime("%Y%m%d_%H%M") + '.txt', 'w') as log:
		log.write('CLAM Re-aligner ' + __version__ + '\n')
		log.write('Args:\n' + '\n'.join(args) + '\n')
		log.write('resume: ' + str(options.resume) + '\n')
		log.write('verbose: ' + str(options.verbose) + '\n')
		log.write('output_dir: ' + str(options.output_dir) + '\n')
		log.write('tmp_dir: ' + str(options.tmp_dir) + '\n')
		log.write('window_size: ' + str(options.window_size) + '\n')
		log.write('max_multihits: ' + str(options.max_multihits) + '\n')
		log.write('is_stranded: ' + str(options.is_stranded) + '\n')
		log.write('max-gap: ' + str(options.max_gaps) + '\n')
		#log.write('gtf: ' + str(options.gtf) + '\n')
		#if len(args)>1:
		#	log.write('cov_site_min: ' + str(options.cov_site_min) + '\n')
		#	log.write('cov_gene_min: ' + str(options.cov_gene_min) + '\n')
	return

def read_rm_out(filename):
	"""
	sub-routine to check nodes with finished EM. Called if resume is turned on.
	"""
	nodes_finished=set()
	with open(filename, 'r') as f:
		for line in f:
			ele=line.split('\t')
			if len(ele) < 7:
				continue
			node_name, in_seen=ele[3].split('|')
			if in_seen=='T':
				nodes_finished.add(node_name)
	return nodes_finished

def read_cufflinks(filename):
	"""
	sub-routine
	"""
	gene_fpkm={}
	with open(filename,'r') as f:
		f.readline()
		for line in f:
			ele=line.split('\t')
			gene_fpkm[ele[0]] = float(ele[9])
	return gene_fpkm
	
def search_node_subg(node, location_to_reads, read_to_locations, seen):
	"""
	Extract the complete independent subgraph given a node, and add all subgraph nodes into the record.
	Returns the nodes in this subgraph and updated record.
	"""
	if node in seen:
		return None, seen
	tmp_net=set()
	queue=[node]
	while(len(queue)>0):		
		map(tmp_net.add, queue)
		map(seen.add, queue)
		new_queue=deque([])
		for x in queue:
			x_reads=location_to_reads[x]
			new_queue.extend([next_node for x_read in x_reads for next_node in read_to_locations[x_read] if not next_node in tmp_net])
		queue=list(set(new_queue))
	subg=list(set(tmp_net))
	return subg, seen

class Bit:
	"""
	Binary Indexed Tree to store values in genomic intervals.
	Implementation modified from http://www.geeksforgeeks.org/binary-indexed-tree-or-fenwick-tree-2/
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


def construct_track_lite(subgraph, location_to_reads, read_to_locations):
	"""
	Construct BIT for each genomic region / node.
	Returns a node-track dictionary and a dictionary for multi-mapped reads.
	"""
	node_track = {}
	total_len = 0
	obs_reads = deque([])
	for node in subgraph:
		chr, strand, start, end = node.split(':')
		start, end = int(start), int(end)
		this_len = end - start + 1
		node_track[node] = Bit(this_len)
		node_reads=location_to_reads[node]
		obs_reads.extend(node_reads)
	
	obs_reads = list(set(obs_reads))
	multi_reads_weights = defaultdict(dict)
	for i in range(len(obs_reads)):
		read = obs_reads[i]
		read_nodes = read_to_locations[read]
		read_score = 1.0 / len(read_nodes)
		is_multi = True if read_score!=1 else False
		for nd in read_nodes:
			chr, strand, nds, nde = nd.split(':')
			nds, nde = int(nds), int(nde)
			rds, rde = [int(x) for x in read_nodes[nd].split('\t')]
			rlen = rde - rds + 1
			read_ind = rds - nds + 1
			read_end = read_ind + rlen - 1 
			read_center = int(np.ceil( (read_ind + read_end) / 2.0))
			node_track[nd].add(read_center, read_score)
			if is_multi:
				multi_reads_weights[read][nd]=[read_score, read_ind, read_end]
	return(node_track, multi_reads_weights)	

def runEM(node_track, multi_reads_weights, w=50, epsilon=1e-6, max_iter=100, verbose=True):
	"""
	EM implementation for re-assigning multi-mapped reads, given the compatibility matrix of a subgraph.
	"""
	iter=0
	residue=1
	while iter < max_iter and residue > epsilon:
		residue = 0
		reweight=defaultdict(dict)
		# calculate re-distribute probability: M-step
		for read in multi_reads_weights:
			for nd in multi_reads_weights[read]:
				sz=node_track[nd].size-1
				old_score, rind, rend = multi_reads_weights[read][nd]
				rcenter = int(np.ceil((rind + rend) / 2.0))
				reweight[read][nd] = max( 0, node_track[nd].sum(min(sz, rcenter + w)) - node_track[nd].sum(max(0,rcenter - w)) )
		# update track; E-step
		for read in reweight:
			dn=sum([reweight[read][x] for x in reweight[read]])
			if dn==0:
				print >> sys.stdout, 'error occured.'
				dn=1
			for nd in reweight[read]:
				old_score, rind, rend = multi_reads_weights[read][nd]
				rcenter = int(np.ceil((rind+rend)/2.0))
				new_score = reweight[read][nd] / float(dn)
				node_track[nd].add(rcenter, new_score - old_score)
				residue += (old_score - new_score)**2
				multi_reads_weights[read][nd][0] = new_score
		if verbose and (not iter % 10 or iter == max_iter):
			print_time_stamp('Iter %d, residue = %f' % (iter, residue))
		iter += 1
	return multi_reads_weights

	
def get_genomic_regions(bamfile, distance, verbose, tmp_dir, is_stranded, min_unique_reads, multiread_set):
	"""
	Construct genomic regions by collapsing reads; also construct read-location dict while iterating the bam file.
	"""
	pileup=defaultdict(list)
	location_to_reads=defaultdict(list)
	read_to_locations=defaultdict(dict)
	head=bamfile.header['SQ']
	chrs=[x['SN'] for x in head]	
	print_time_stamp('Finding genomic regions with more than ' + str(min_unique_reads) + ' unique reads.')
	for chr in chrs:
		chr_aln=[x for x in bamfile.fetch(chr)]
		tmp_cluster_pos=[0, 0, 0]
		tmp_cluster_neg=[0, 0, 0]
		tags_pos=[]
		tags_neg=[]
		if verbose:
			print_time_stamp('finding genomic regions on ' + chr + '..')
		for read in chr_aln:
			read_len = read.qlen if read.qlen > 0 else read.positions[-1] - read.positions[0] + 1
			if read_len > 100 or read.positions[-1] - read.positions[0] > 100:  # possibly a junction read, discard
				continue
			pos=[int(x) for x in read.positions]
			if read.is_reverse and is_stranded:  # add here to avoid negative strand if not stranded library
				if (pos[0] + pos[-1])/2 - tmp_cluster_neg[2] <= distance:  # compute distance using read centers.
					tmp_cluster_neg[1] = max(pos[-1], tmp_cluster_neg[1])
					tmp_cluster_neg[2] = (pos[0] + pos[-1])/2
					tags_neg.append([read.qname, read.positions[0], read.positions[-1]])
				else:
					if len(tags_neg) > 1 and sum([1 for x in tags_neg if not x[0] in multiread_set]) >= min_unique_reads:
						node_name = chr + ':-:' + str(tmp_cluster_neg[0]) + ':' + str(tmp_cluster_neg[1])
						pileup[chr + ':-'].append([tmp_cluster_neg[0], tmp_cluster_neg[1]])
						location_to_reads[node_name].extend([x[0] for x in tags_neg])
						for x_qname, x_pos0, x_pos1 in tags_neg:
							read_to_locations[x_qname].update({node_name : str(x_pos0) + '\t' + str(x_pos1) })
					tmp_cluster_neg=[pos[0], pos[-1], (pos[0] + pos[-1])/2]
					tags_neg=[ [read.qname, read.positions[0], read.positions[-1]] ]
			else:
				if (pos[0] + pos[-1])/2 - tmp_cluster_pos[2] <= distance:
					tmp_cluster_pos[1]=max(pos[-1], tmp_cluster_pos[1])
					tmp_cluster_pos[2] = (pos[0] + pos[-1])/2
					tags_pos.append([read.qname, read.positions[0], read.positions[-1]])
				else:
					if len(tags_pos) > 1 and sum([1 for x in tags_pos if not x[0] in multiread_set]) >= min_unique_reads:
						node_name = chr + ':+:' + str(tmp_cluster_pos[0]) + ':' + str(tmp_cluster_pos[1])
						pileup[chr + ':+'].append([tmp_cluster_pos[0], tmp_cluster_pos[1]])
						location_to_reads[node_name].extend([x[0] for x in tags_pos])
						for x_qname, x_pos0, x_pos1 in tags_pos:
							read_to_locations[x_qname].update({node_name : str(x_pos0) + '\t' + str(x_pos1) })
					tmp_cluster_pos=[pos[0], pos[-1], (pos[0] + pos[-1])/2]
					tags_pos=[ [read.qname, read.positions[0], read.positions[-1]] ]		
	print_time_stamp('Save genomic regions to file.')
	pickle.dump(pileup, open(tmp_dir + '/genomic_regions_%s_%s.pdata' % (str(distance), str(min_unique_reads)),'wb'), -1)
	pickle.dump(location_to_reads, open(tmp_dir + '/loc2read.pdata','wb'), -1)
	pickle.dump(read_to_locations, open(tmp_dir + '/read2loc.pdata','wb'), -1)
	return pileup, location_to_reads, read_to_locations

def update_by_filter(genomic_regions, read_to_locations, location_to_reads, preserved_nodes):
	"""
	Remove nodes that don't pass the filtering.
	Returns filtered nodes/genomic regions.
	"""
	new_genomic_regions = defaultdict(list)
	new_read_to_locations = defaultdict(dict)
	new_location_to_reads = defaultdict(list)
	for chr_strand in genomic_regions:
		for start, end in genomic_regions[chr_strand]:
			node = chr_strand + ':' + str(start) + ':' + str(end)
			if node in preserved_nodes:
				new_genomic_regions[chr_strand].append([start, end])
	for read in read_to_locations:
		for node in read_to_locations[read]:
			if node in preserved_nodes:
				new_read_to_locations[read].update({node:read_to_locations[read][node]})
	for location in location_to_reads:
		if location in preserved_nodes:
			new_location_to_reads[location]=location_to_reads[location]
	return(new_genomic_regions, new_read_to_locations, new_location_to_reads)
	
def filter_by_covr(cov_filename, cov_site_min, cov_gene_min, gtffile, genomic_regions, tmp_dir):
	"""
	Wrapper for filtering nodes. Filter based on minimum unique read counts and minimum gene expression.
	"""
	node_list = [ chr_strand+':'+str(start)+':'+str(end) for chr_strand in genomic_regions for start,end in genomic_regions[chr_strand] ]
	covfile=pysam.Samfile(cov_filename)
	gene_preserved = deque()
	site_preserved = deque()
	
	if cov_site_min > 0:
		k = 0
		for chr_strand in genomic_regions:
			chr, strand = chr_strand.split(':')
			print_time_stamp('filtering site: '+chr_strand)
			for start, end in genomic_regions[chr_strand]:
				k += 1
				if not k % 10000:
					print_time_stamp('filtering site count: ' + str(k) + '/' + str(len(node_list)))
				node = chr_strand + ':' + str(start) + ':' + str(end)
				num_reads = sum([1 for x in covfile.fetch(chr, int(start), int(end)) if x.pos>start and x.pos<end])
				if num_reads >= cov_site_min:
					site_preserved.add(node)
	else:
		site_preserved=set(node_list)
	
	if cov_gene_min>0:
		genomic_regions_list = [ (chr_strand.split(':')[0], int(start), int(end), chr_strand+':'+':'.join([str(start),str(end)]), 'X', chr_strand.split(':')[1] ) for chr_strand in genomic_regions for start,end in genomic_regions[chr_strand]]
		genomic_regions_bed = pybedtools.BedTool(genomic_regions_list)
		gtf = pybedtools.BedTool(gtffile)
		overlap_transcripts = genomic_regions_bed.intersect(gtf, wo=True, s=True)
		overlap_transcripts.saveas(tmp_dir + '/genomic_regions.gtf.bed')
		total=len(overlap_transcripts)
		pybedtools.cleanup()
		del overlap_transcripts
		del gtf
		del genomic_regions_list
		
		cov_scale=sum([int(x.split('\t')[2]) for x in pysam.idxstats(cov_filename).split('\n') if len(x)>0]) / 1000000.0
		
		#gene_fpkm=read_cufflinks('/u/home/f/frankwoe/scratch/Ule_RNAseq_hnRNPC/cufflinks_output_star/genes.fpkm_tracking')
		gene_rpkm={}
		k = 0
		f = open(tmp_dir + '/genomic_regions.gtf.bed','r')
		for ele in f:
			line = ele.split()
			k += 1
			if not k % 10000:
				print_time_stamp('filtering gene RPKM: ' + str(k) + '/' + str(total))
			node=line[3]
			#if not node in site_preserved:
			#	continue
			gene_id = line[9]
			#RPKM = gene_fpkm[gene_id] if gene_id in gene_fpkm else 0
			if gene_id in gene_rpkm:
				RPKM = gene_rpkm[gene_id]
			else:
				chr, start, end = line[6], line[7], line[8]
				transcript_count = covfile.count(chr, int(start), int(end))
				block_sizes = [int(x) for x in line[16].split(',') if x!='']
				gene_len = sum(block_sizes) / 1000.0
				RPKM = transcript_count / cov_scale / gene_len
				gene_rpkm[gene_id]=RPKM
			
			if RPKM >= cov_gene_min:
				gene_preserved.append(node)
		gene_preserved=set(gene_preserved)
		f.close()
	else:
		gene_preserved=set(node_list)
	
	return gene_preserved.intersection(site_preserved)
	

def make_dependency_dict(tmp_dir, is_stranded):
	sort_cmd = 'LC_COLLATE=C sort -k1,1 -k2,2n %s/reg.bed > %s/reg_sorted.bed' % (tmp_dir, tmp_dir)
	subprocess.call(sort_cmd, shell=True)
	strand_opt = '-s' if is_stranded else ''
	multi_reg_cmd = 'bedtools intersect %s -a %s/reg_sorted.bed -b %s/multi.bed -wo -sorted > %s/reg_multi.bed' % (strand_opt, tmp_dir, tmp_dir, tmp_dir)
	subprocess.call(multi_reg_cmd, shell=True)
	uni_reg_cmd = 'bedtools intersect %s -a %s/reg_sorted.bed -b %s/unique.bed -wo -sorted > %s/reg_unique.bed' % (strand_opt, tmp_dir, tmp_dir, tmp_dir)
	subprocess.call(uni_reg_cmd, shell=True)
	
	
def print_time_stamp(msg):
	"""
	Reporter function for logging.
	"""
	current_time='[' + strftime("%Y-%m-%d %H:%M:%S") + '] '
	print >> sys.stderr, current_time + msg

def filter_multihits(filename, max_hits, verbose, tmp_dir):
	"""
	Pre-processing function for cleaning up the input bam file.
	"""
	if verbose:
		print_time_stamp('filtering multi-mapped up to %d hits.' % max_hits)
	multiread_set=set()
	subprocess.call("samtools view -h %s | awk '{if($2 !~ /_/ && $3 !~ /_/) {print}}' | samtools view -bS - > %s/filter_random.bam" % (filename, tmp_dir), shell=True)
	oldfile=pysam.Samfile(tmp_dir + '/filter_random.bam','rb')
	new_filename=os.path.abspath(tmp_dir + '/filter%d.bam' % max_hits)
	sorted_filename=os.path.abspath(tmp_dir + '/filter%d.sorted.bam' % max_hits)
	newfile=pysam.Samfile(new_filename, 'wb', template=oldfile)
	for aligned_read in oldfile:
		try:
			if aligned_read.opt("NH") < max_hits:
				newfile.write(aligned_read)
				if aligned_read.opt("NH")>1:
					multiread_set.add(aligned_read.qname)
		except KeyError:
			newfile.write(aligned_read)
	oldfile.close()
	newfile.close()
	sort_cmd='samtools sort -T %s/ -o %s %s' % (tmp_dir, sorted_filename, new_filename)
	index_cmd='samtools index %s' % sorted_filename
	subprocess.call(sort_cmd, shell=True)
	subprocess.call(index_cmd, shell=True)
	subprocess.call('rm %s/filter_random.bam %s' % (tmp_dir, new_filename), shell=True)
	pickle.dump(multiread_set, open(tmp_dir + '/multiread_set.pdata', 'wb'), -1)
	return(sorted_filename, multiread_set)

def make_write_content(multi_reads_weights):
	"""
	Make the content for output files based on results in memory.
	"""
	content=''
	for read in multi_reads_weights:
		for node in multi_reads_weights[read]:
			score, rds, rde = multi_reads_weights[read][node]
			chr, strand, nds, nde = node.split(':')
			read_start = int(nds) + rds - 1
			read_len = int(rde) - int(rds) + 1
			read_end = read_start + read_len - 1
			content += '\t'.join([chr, str(read_start), str(read_end), read, str(score), strand]) + '\n'
	return content

if __name__=='__main__':
	main()
