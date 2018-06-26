#!/usr/bin/env python

"""This peak-caller script is part of the CLAM pipeline.

This subcommand will call peaks using permutation by randomly placing reads along the gene.
More details about the permutation procedure is described in our NAR paper.

Example run:
	```
	CLAM permutation_callpeak -i path/to/outdir/unique.sorted.bam path/to/outdir/realigned.sorted.bam \
	-o path/to/peaks/outdir -p 8 \
	--gtf path/to/gencode.v19.annotation.gtf
	```
Author:
	Zijun Zhang <zj.z@ucla.edu>

Tested under python 2.7.3
"""

__author__ = 'Zijun Zhang'
__version__ = '1.1.3'
__email__ = 'zj.z@ucla.edu'


import os
import sys
from collections import defaultdict
from statsmodels.sandbox.stats.multicomp import multipletests
import logging
import bisect
import random
import pysam
import re
from multiprocessing import Pool


###setup logger
logger = logging.getLogger('CLAM.permutation_peakcaller')
###

def parser(args):
	"""The main wrapper for CLAM peak-caller.
	"""
	# logging info
	logger.info('start')
	logger.info('run info: %s'%(' '.join(sys.argv)))
	# some back-reference (to v1.0.0) parameters here
	random_state = args.random_state
	merge_size = args.merge_size
	output_dir = os.path.abspath(args.out_dir)
	nthread = args.nthread
	max_iter = 200
	
	# read in gtf gene annotations
	gene_annot = read_gtf(args.gtf_fp)
	
	# read in GTF
	gene_list = gene_annot.keys()
	child_gene_list = [x for x in chunkify(gene_list, nthread)]	
	
	# call peaks
	unibam_file=args.in_bam[0]
	multibam_file=args.in_bam[1]	
	if nthread>1:
		pool = Pool(processes=args.nthread)
		assert len(args.in_bam)==2
		tid_to_qval_compact = pool.map(
			_child_get_permutation_fdr, 
			[ (unibam_file, multibam_file, child_gene_list[i], gene_annot, args.qval_cutoff, max_iter, ~args.unstranded, 'fdr', random_state)
				for i in range(args.nthread) 
			])

		pool.terminate()
		pool.join()

		unique_tid_to_qval, combined_tid_to_qval = unpack_tid_to_qval(tid_to_qval_compact)
	else:
		unique_tid_to_qval, combined_tid_to_qval = _child_get_permutation_fdr(
				(unibam_file, multibam_file, gene_list, gene_annot, args.qval_cutoff, max_iter, ~args.unstranded, 'fdr', random_state)
			)
	
	#pickle.dump(unique_tid_to_qval, open(tmp_dir+'/unique_to_qval.pdata','wb'), -1)
	#pickle.dump(combined_tid_to_qval, open(tmp_dir+'/combined_to_qval.pdata','wb'), -1)
	merge_peaks = merge_peaks_singleNucl
	#if args.merge_method==1:
	#	merge_peaks=merge_peaks_singleNucl
	#	mm='singleNucl'
	#elif args.merge_method==2:
	#	merge_peaks=merge_peaks_broadPeak
	#	mm='broadPeak'
	#else:
	#	merge_peaks=merge_peaks_singleNucl
	#	mm='unknown selection, using default singleNucl'
		
	
	unique_peaks=merge_peaks(unique_tid_to_qval, merge_size, args.qval_cutoff)
	combined_peaks=merge_peaks(combined_tid_to_qval, merge_size, args.qval_cutoff)
	
	# write peak-calling results to file.
	narrowPeak_formatter = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t.\t%.3e\t.\n"
	## chr start end name unique/combined strand signalValue pVal qVal peak
	with open(output_dir + '/all_permutation_peaks.bed', 'w') as f:
		for peak in unique_peaks:  # peak = ['chr\tstart\tend\tstrand', 'height\tqval\t', tid]
			if args.extend is None:
				wt_loc=peak[0]
			else:
				wt_loc=extend_peak_region(peak[0], args.extend)
			#f.write(wt_loc + '\t' + '\t'.join([str(x) for x in peak[1]]) + '\t' + peak[2] + '\tunique\n')
			chr, start, end, strand = wt_loc.split('\t')
			_, signal_qval, gene_name = peak
			signal, qval = signal_qval
			f.write( narrowPeak_formatter % (chr, start, end, gene_name, 'unique', strand, signal, qval) )
		for peak in combined_peaks:
			if args.extend is None:
				wt_loc=peak[0]
			else:
				wt_loc=extend_peak_region(peak[0], args.extend)
			#f.write(wt_loc + '\t' + '\t'.join([str(x) for x in peak[1]]) + '\t' + peak[2] + '\tcombined\n')
			chr, start, end, strand = wt_loc.split('\t')
			_, signal_qval, gene_name = peak
			signal, qval = signal_qval
			f.write( narrowPeak_formatter % (chr, start, end, gene_name, 'combined', strand, signal, qval) )
	if args.unstranded:
		cmd = ''' sort -k1,1 -k2,2n %s/all_permutation_peaks.bed |awk '{OFS="\t"; print $1,$2,$3,$4":"$7":"$9,$5,$6}'| bedtools merge -d -1 -i stdin -c 4,5,6 -o collapse,collapse,distinct  > %s''' % (output_dir, os.path.join(output_dir,'narrow_peak.permutation.bed') )
	else:
		cmd = ''' sort -k1,1 -k2,2n %s/all_permutation_peaks.bed |awk '{OFS="\t"; print $1,$2,$3,$4":"$7":"$9,$5,$6}'| bedtools merge -s -d -1 -i stdin -c 4,5,6 -o collapse,collapse,distinct  > %s''' % (output_dir, os.path.join(output_dir,'narrow_peak.permutation.bed') )
	os.system( cmd )
	logger.info('end')
	return


def chunkify(a, n):
	"""
	Separate a list (a) into consecutive n chunks.
	Returns the chunkified index
	"""
	k, m = len(a) / n, len(a) % n
	return (a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in xrange(n))

	
def unpack_tid_to_qval(compact):
	"""
	Unpacks the returned values from multi-processing.
	"""
	unique_tid_to_qval=defaultdict(list)
	combined_tid_to_qval=defaultdict(list)
	for item in compact:
		unique, combined = item[0], item[1]
		for tid in combined:
			if len(unique[tid])>0:
				unique_tid_to_qval[tid]=unique[tid]
			if len(combined[tid])>1:
				combined_tid_to_qval[tid]=combined[tid]
	return unique_tid_to_qval,combined_tid_to_qval

	
def _child_get_permutation_fdr(
		(unibam_file, multibam_file, child_gene_list, gene_annot, 
		pval_cutoff, max_iter, is_stranded, correction_method,seed)
	):
	"""
	General permutation wrapper for a list of genes. Gets called by multi-processing generated by Pool()
	Returns packed FDRs from each child process.
	"""
	random.seed(seed)
	
	unique_tid_to_qval=defaultdict(list)
	combined_tid_to_qval=defaultdict(list)
	
	unibam=pysam.Samfile(unibam_file, 'rb')
	multibam=pysam.Samfile(multibam_file, 'rb')
	
	pid = os.getpid()
	tot = len(child_gene_list)
	
	for i in range(len(child_gene_list)):
		if not i % 200:
			logger.debug('pid %s : %i / %i (%.2f%%)'% (pid, i, tot, float(i)/float(tot)*100))
		gene_name = child_gene_list[i]
		gene = gene_annot[gene_name]
		chr, start, end, strand, tid = gene[0:5]
		unique_reads = read_tid_frag_from_bam(gene, unibam, is_stranded, True)
		multi_reads = read_tid_frag_from_bam(gene, multibam, is_stranded, False)
		
		this_unique_to_qval = do_permutation(gene, unique_reads, max_iter, pval_cutoff, correction_method)
		this_combined_to_qval = do_permutation(gene, unique_reads+multi_reads, max_iter, pval_cutoff, correction_method)
		
		unique_tid_to_qval[tid].extend(this_unique_to_qval)
		combined_tid_to_qval[tid].extend(this_combined_to_qval)
	unibam.close()
	multibam.close()
	return unique_tid_to_qval, combined_tid_to_qval


def do_permutation(transcr, read_transcript, max_iter, pval_cutoff, correction_method):	
	"""
	Permutes the reads along a given gene length, sub-routine that get called by get_permutation_fdr(..).
	Returns the locally corrected p-values for each observed height on the given gene.
	"""
	chr, tstart, tend, strand, tid = transcr[0:5]
	tid_length=tend-tstart+1
	obs_heights_count=count_pileup_heights(tid_length, read_transcript)
	
	tid_to_qval=[]
	
	rand_heights_dist=defaultdict(int)
	rand_sum=0
	# need to account for the 'observed' data, since permutation tests should never report p-value as 0. 3/22/16
	for i in obs_heights_count:
		if i==0:
			continue
		else:
			rand_heights_dist[int(i)]+=1
			rand_sum+=1
	for B in range(max_iter):
		new_heights_count=permutate_heights(tid_length, read_transcript)
		for i in new_heights_count:
			if i==0:
				continue
			else:
				rand_heights_dist[i]+=1
				rand_sum+=1
	height_to_pval={}
	for h in set(obs_heights_count):
		if h < 1:
			continue
		else:
			lefter=0
			for j in range(int(h), max(rand_heights_dist)+1):
				lefter+=rand_heights_dist[j]
			height_to_pval[h]=lefter/float(rand_sum)
	pval_list=[]
	for i in obs_heights_count:
		if i<1:
			continue
		pval_list.append(height_to_pval[i])
	if len(pval_list)<=1:
		return []
	
	qval_list=multipletests(pval_list, method='fdr_bh')[1]
	#if correction_method==2 or correction_method.lower()=='fdr':
	#	qval_list=multipletests(pval_list, method='fdr_bh')[1]
	#else:
	#	qval_list=[min(x*(len(set([int(y) for y in height_to_pval if y!=0]))), 1.0) for x in pval_list]
	
	ind=0
	last_height=0
	for j in range(len(obs_heights_count)):
		this_height=obs_heights_count[j]
		if this_height<1:
			last_height=0
			continue
		if qval_list[ind] <= pval_cutoff:
			if this_height==last_height:
				chr, last_start, last_end, last_strand, last_height, last_qval=tid_to_qval[-1]
				tid_to_qval[-1]=[chr, last_start, tstart+j+1, strand, last_height, last_qval]
			else:
				tid_to_qval.append([chr, tstart+j, tstart+j+1, strand, obs_heights_count[j], qval_list[ind]])  # chr, start, end, strand, height, this_qval
				last_height=this_height
		ind+=1
	return tid_to_qval


def heights_to_dist(rand_heights):
	"""
	sub-routine 
	"""
	rand_heights_dist=defaultdict(int)
	rand_sum=0
	for new_heights_count in rand_heights:
		for i in new_heights_count:
			if i==0:
				continue
			else:
				rand_heights_dist[i]+=1
				rand_sum+=1
	return rand_heights_dist, rand_sum


def permutate_heights(tlen, reads):
	"""
	Sub-routine for do_permutation(...)
	Randomly allocate the read locations.
	"""
	loc_heights=[0] * tlen
	for id, pos, read_len, score in reads:
		if score<1 and random.random() > score:
			continue
		rand_pos=random.randint(1, max(1, tlen-read_len))
		for i in range(rand_pos, min(rand_pos + read_len, tlen)):
			loc_heights[i]+=1
	return loc_heights


def count_pileup_heights(tlen, reads):
	"""
	Sub-routine for do_permutation(...)
	Counts the distribution of pile-up heights for a given gene/permutation
	"""
	loc_heights=[0] * tlen
	for id, pos, read_len, score in reads:
		for i in range(pos, min(pos+read_len-1, tlen)):
			loc_heights[i]+=score
	return loc_heights


def merge_peaks_broadPeak(transcript_to_qval, merge_size, pval_cutoff):
	"""
	Merge called peaks on a gene using option 2, 
	i.e. if two peaks close to each other, region
	between two peaks are also called as peaks
	Retuns a list of merged peaks.
	"""
	peaks=[]
	last_qval=[0,1]
	for tid in transcript_to_qval:
		init=True
		for chr, start, end, strand, height, this_qval in transcript_to_qval[tid]:
			loc=[chr, str(start), str(end), strand]
			this_qval=[height, this_qval]  # this_qval=[height, qval] so that when qval=0, we can compare height
			if  this_qval[1] > pval_cutoff:
				continue
			if init:
				last_qval=this_qval
				last_pos=[start, end]
				last_loc=loc
				last_chr=chr
				write_out=False
				init=False
				continue
			if int(start) - int(last_pos[1]) > merge_size:
				write_out=True
			else:
				last_pos=[last_pos[0], end]
				last_qval=this_qval if last_qval[0]<this_qval[0] else last_qval
				last_loc[2]=str(end)
				write_out=False
				
			if write_out and last_qval[1] < pval_cutoff:
				peaks.append(['\t'.join(last_loc), last_qval, tid])
				last_qval=this_qval
				last_pos=[start, end]
				last_loc=loc
				last_chr=[chr, str(start), str(end), strand]
				write_out=False
		if last_qval[1] < pval_cutoff:
			peaks.append(['\t'.join(last_loc), last_qval, tid])
	return peaks


def read_gtf(fn):
	"""read in the gene annotation from GTF file
	"""
	logger.info('read GTF from "%s" '% fn)
	gene_annot = {}
	with open(fn, 'r') as f:
		for line in f:
			if line.startswith('#'):
				continue
			ele = line.strip().split('\t')
			if ele[2] != 'gene':
				continue
			chr, start, end, strand = ele[0], int(ele[3]), int(ele[4]), ele[6]
			try:
				gene_id = re.search(r'gene_id "(.+?)"', ele[-1]).group(1)
			except AttributeError:
				continue
			gene_annot[gene_id] = [chr, start, end, strand, gene_id]
	return gene_annot


def merge_peaks_singleNucl(transcript_to_qval, merge_size, pval_cutoff):
	"""
	Merge called peaks on a gene using option 1 
	(default), i.e. if two peaks close to each other, 
	only pick the most significant one peak
	Retuns a list of merged peaks.
	"""
	peaks=[]
	last_qval=[0,1]
	for tid in transcript_to_qval:
		init=True
		for chr, start, end, strand, height, this_qval in transcript_to_qval[tid]:
			loc='\t'.join([chr, str(start), str(end), strand])
			this_qval=[height, this_qval]  # this_qval=[height, qval] so that when qval=0, we can compare height
			if  this_qval[1] > pval_cutoff:
				continue
			if init:
				last_qval=this_qval
				last_pos=[start, end]
				last_loc=loc
				last_chr=chr
				write_out=False
				init=False
				continue
			if last_chr == chr:
				if abs( int(start) - int(last_pos[0]) ) > merge_size:
					write_out=True
				elif last_qval[0] < this_qval[0]:
					last_pos=[start, end]
					last_qval=this_qval
					last_loc=loc
					write_out=False
			else:
				write_out=True
				
			if write_out and last_qval[1] < pval_cutoff:
				#peaks[last_loc]=last_qval
				peaks.append([last_loc, last_qval, tid])
				last_qval=this_qval
				last_pos=[start, end]
				last_loc=loc
				last_chr=chr
				write_out=False
		if last_qval[1] < pval_cutoff:
			peaks.append([last_loc, last_qval, tid])
	return peaks


def extend_peak_region(loc, target_len):
	"""
	Extends peak symmetrically if peak is smaller than target_len.
	"""
	chr, start, end, strand = loc.split('\t')
	start = int(start)
	end = int(end)
	old_len = end - start
	if old_len > target_len:
		return loc
	else:
		center = int((start + end)/2)
		start = center - int(target_len /2)
		end = center + int(target_len/2)
		return '\t'.join([chr, str(start), str(end), strand])


def read_tid_frag_from_bam(tid, bamfile, is_stranded, is_unique):
	"""
	Use pysam to fetch reads info for a given gene and its loci.
	Returns reads, read weights and its mapped loci.
	"""
	tid_reads=[]
	#gene, chr, strand, start, end=tid
	chr, start, end, strand, gene = tid[0:5]
	if strand=='-':
		is_reverse=True
	else:
		is_reverse=False
	reads=[x for x in bamfile.fetch(chr, int(start), int(end)) if x.is_reverse==is_reverse or not is_stranded]
	reads=[x for x in reads if x.pos>=int(start) and x.pos<=int(end)]
	for read in reads:
		if is_unique:
			try:
				opt_NH=read.opt('NH')
				if opt_NH > 1:
					continue
			except:
				pass
			score=1
		else:
			try:
				opt_AS=read.opt('AS')
				if isinstance(opt_AS, float):
					score=opt_AS
				else:
					continue
			except:
				continue
		#try:
		#	read_length = read.qlen
		#except:
		#	read_length = read.positions[-1] - read.positions[0] + 1
		read_length = read.positions[-1] - read.positions[0] + 1
		
		if (not 'N' in read.cigarstring) and \
			(read.pos-start>=0 and read_length<500):    # to avoid junction reads
				tid_reads.append([read.qname, read.pos-start, read_length, score])
	return tid_reads


