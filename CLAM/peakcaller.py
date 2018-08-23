#!/usr/bin/env python

"""This peak-caller script is part of the CLAM pipeline.

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
import pysam
import logging
import numpy as np
from collections import defaultdict
import re
from scipy.stats import fisher_exact, poisson, chi2
import scipy.optimize as optimize
from tqdm import tqdm
import datetime
from .stats import ztnb_em, bin_test_alternatives
from multiprocessing import Pool
import argparse as ap
import inspect


### get logger
###
logger = logging.getLogger('CLAM.Peakcaller')
###

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
				#gene_id = re.search(r'gene_id "(.+?)"', ele[-1]).group(1)
				## a new regular expression generalizable to 
				## universal situations
				gene_id = re.findall(r"(\w+)", ele[-1])[1]
			except AttributeError:
				continue
			gene_annot[gene_id] = [chr, start, end, strand, gene_id]
	return gene_annot


def count_gene_read_tags(bam_list, gene, is_unique=True, unstranded=False):
	""" count the tagger positions for all reads in a given genomic interval
	Args:
	Returns:
	"""
	chr, start, end, strand = gene[0:4]
	# placeholder for interval: 'num of replicate' x 'interval length'
	interval = np.zeros( (len(bam_list), end-start+1) )
	is_reverse = True if strand=='-' else False
	# construct the (tag, score) pairs
	for i in range(len(bam_list)):
		bam = bam_list[i]
		if is_unique:
			read_tags = [ (x.opt('RT'), 1.0) for x in bam.fetch(chr, start, end) \
				if unstranded or x.is_reverse==is_reverse]
		else:
			#read_tags = [ (x.opt('RT'), x.opt('AS')) for x in bam.fetch(chr, start, end) \
			#	if unstranded or x.is_reverse==is_reverse]
			## UPDATE 1.22.2018: I noted in certain genes with rRNA, when considering 
			## multi-mapped reads, the memory usage will significantly increase just for
			## that single gene, rendering `MemoryError` even on single threads
			## This is especially a problem for RIP-seq controls.
			## For now, let's set a very high threshold for maximum number of tags in any
			## given gene. This should not cause any problem in CLIP-seq, and should
			## not affect too much results.
			## If this happens, will log a INFO-level message in the log
			counter = 0
			read_tags = []
			for x in bam.fetch(chr, start, end):
				counter += 1
				if unstranded or x.is_reverse==is_reverse:
					read_tags.append([x.opt('RT'), x.opt('AS')])
				if counter > 0.5*10**8:
					logger.info("too many mreads (>50million) on %s, probably rRNA gene - discarded"%(':'.join([chr,strand,str(start),str(end)])))
					return interval
		for tag in read_tags:
			if tag[0]<start or tag[0]>=end:
				continue
			interval[i, tag[0]-start] += tag[1]
	return interval


def	bin_interval_counts(interval, binsize=50):
	bins = np.zeros( ( interval.shape[0], int(np.ceil(interval.shape[1]/float(binsize))) ) )
	for i in range(bins.shape[1]):
		for j in range(interval.shape[0]):
			start, end = i*binsize, (i+1)*binsize
			bins[j, i] = np.sum(interval[j, start:end])
	return bins


def test_bin_negbinom(intv_bin_ip, intv_bin_con, with_control=True, correction_method='fdr_bh'):
	"""DOCSTRING
	Args
	Returns
	"""
	def _par_to_vec(par, data, is_constrained):
		if is_constrained:
			beta = par[0]
			mu_vec = par[1::]
			delta = 0
		else:
			beta, delta = par[0], par[1]
			mu_vec = par[2::]
		ip_counter = data['this_ip'].shape[0]
		con_counter = data['this_con'].shape[0]
		mu0 = np.asarray(mu_vec[0:con_counter])
		mu1 = np.asarray(mu_vec[con_counter::])
		lamb1_this = np.exp(mu1 + beta + delta)
		lamb1_others = np.exp(mu1)
		lamb0_this = np.exp(mu0 + beta)
		lamb0_others = np.exp(mu0)
		return (lamb1_this, lamb1_others, lamb0_this, lamb0_others)
	
	def _negative_binom_logpmf(y, mu, alpha):
		y  = np.asarray(y)
		ll = np.empty(len(y))
		for i in range(len(y)):
			alpha_inv = 1.0/alpha[i]
			alpha_mu = float(alpha[i] * mu[i])
			ll[i] = y[i]* np.log(alpha_mu/(1+alpha_mu))- \
				alpha_inv*np.log(1+alpha_mu)
		return ll
		
	def _neg_loglik_unconstrain(par, data):
		(l1, l2, l3, l4) = _par_to_vec(par, data, False)
		ll = np.sum( _negative_binom_logpmf(data['this_ip'], mu=l1, alpha=alpha_ip_vec))
		ll += np.sum( _negative_binom_logpmf(data['others_ip'], mu=l2, alpha=alpha_ip_vec))
		ll += np.sum( _negative_binom_logpmf(data['this_con'], mu=l3, alpha=alpha_con_vec))
		ll += np.sum( _negative_binom_logpmf(data['others_con'], mu=l4, alpha=alpha_con_vec))
		return -ll
	
	def _neg_loglik_constrain(par, data):
		(l1, l2, l3, l4) = _par_to_vec(par, data, True)
		ll = np.sum(_negative_binom_logpmf(data['this_ip'], mu=l1, alpha=alpha_ip_vec)) + \
			np.sum(_negative_binom_logpmf(data['others_ip'], mu=l2, alpha=alpha_ip_vec)) + \
			np.sum(_negative_binom_logpmf(data['this_con'], mu=l3, alpha=alpha_con_vec)) + \
			np.sum(_negative_binom_logpmf(data['others_con'], mu=l4, alpha=alpha_con_vec))
		return -ll
	
	# initialize placeholders
	intv_counter = intv_bin_ip.shape[1]
	assert intv_counter == intv_bin_con.shape[1]
	binscore = np.empty(intv_counter)
	binsignal = np.empty(intv_counter, dtype='S50')
	alpha_ip_vec = np.empty(intv_bin_ip.shape[0])
	alpha_con_vec = np.empty(intv_bin_con.shape[0])
	ip_sum = np.apply_along_axis(np.sum, 1, intv_bin_ip)
	con_sum = np.apply_along_axis(np.sum, 1, intv_bin_con)
	if any(ip_sum==0) or any(con_sum==0):
		return None, None, None
	
	# compute the dispersion parameters
	min_alpha = 0.001 ## small alpha reduces to Poisson
	max_alpha = 0.01   ## large alpha costs loss of power
	if with_control:
		for i in range(intv_bin_con.shape[0]):
			height = ztnb_em.collapse_data(np.floor(intv_bin_con[i,]))
			height[0] = 0
			try:
				ll, mu, alpha = ztnb_em.EM_estim_params(height, max_iter=10, verbose=False)
			except:   ## means truncated loglik=0, not enough data
				alpha = max_alpha
			alpha = max_alpha if alpha>max_alpha else alpha
			alpha = min_alpha if alpha<min_alpha else alpha
			alpha_con_vec[i] = alpha
		## update alpha bounds based on control
		## not until we have a better alpha estimate...  zijun 9.7.2017
		#max_alpha = np.min([ max_alpha, np.max(alpha_con_vec) ] )
		#min_alpha = np.max([ min_alpha, np.min(alpha_con_vec) ] )
		
	
	for i in range(intv_bin_ip.shape[0]):
		height = ztnb_em.collapse_data(np.floor(intv_bin_ip[i,]))
		height[0] = 0
		if np.sum(height.values())>0:
			try:
				ll, mu, alpha = ztnb_em.EM_estim_params(height, max_iter=10, verbose=False)
			except:
				alpha = max_alpha
		else:
			alpha = max_alpha
		## bound IP alpha by background max alpha
		## because IP can be very sparse compared to control,
		## as a result, alpha estimate can be quite off
		alpha = max_alpha if alpha>max_alpha else alpha
		alpha = min_alpha if alpha<min_alpha else alpha
		alpha_ip_vec[i] = alpha
	# use the same dispersion if using *fake* control
	if not with_control:
		alpha_con_vec = np.asarray(alpha_ip_vec)
	
	# this is important:
	# we need to bound the beta parameter, so if the current bin is zero it will be
	# penalized by likelihood
	pseudo_count = 1.
	beta_bound = np.log(pseudo_count/np.min(np.concatenate([ip_sum, con_sum])))
	# perform test on each bin
	for i in range(intv_counter):
		this_ip = intv_bin_ip[:, i]
		others_ip = ip_sum - this_ip
		this_con = intv_bin_con[:, i]
		others_con = con_sum - this_con
		if np.sum(this_ip) == 0:
			binsignal[i], binscore[i] = 'nan', np.nan
			continue
		## to deal with fractional counts, we took the 
		## conservative approach using the floor
		data = {
				'this_ip':np.floor(this_ip),
				'others_ip':np.floor(others_ip),
				'this_con':np.floor(this_con),
				'others_con':np.floor(others_con)
			}
		## constrained likelihood
		res_constrain = optimize.minimize(
				x0=np.ones(1+this_ip.shape[0]+this_con.shape[0]), # beta + mu_vec 
				fun=_neg_loglik_constrain,
				args=(data),
				method='l-bfgs-b',
				bounds = [(beta_bound, abs(beta_bound) )]+[(0.1,10)]*(this_ip.shape[0]+this_con.shape[0]),
				options={'disp':False}
			)
		## unconstrained likelihood
		res_unconstrain = optimize.minimize(
				x0=np.ones(2+this_ip.shape[0]+this_con.shape[0]), # beta + delta + mu_vec
				fun=_neg_loglik_unconstrain,
				args=(data),
				method='l-bfgs-b',
				bounds = [(beta_bound, abs(beta_bound))] + [(-100,100)] + [(0.1,10)]*(this_ip.shape[0]+this_con.shape[0]),
				options={'disp':False}
			)
		
		delta_mle = res_unconstrain.x[1]
		pval = 1 - chi2.cdf(2*(res_constrain.fun - res_unconstrain.fun), 1)
		if delta_mle<0: pval = 1.0
		binscore[i] = pval
		if with_control:
			binsignal[i] = '%.2f'%delta_mle
		else:
			binsignal[i] = ','.join([str(int(x)) if x==int(x) else str(x) for x in this_ip])
	
	# correcting for multiple-testing
	adj = multipletests(binscore[~ np.isnan(binscore)], alpha=0.05, method=correction_method)
	binscore_adj = np.copy(binscore)
	binscore_adj[ ~ np.isnan(binscore) ] = adj[1]
	return binsignal, binscore_adj, binscore


def call_gene_peak(bam_dict, gene, unique_only=False, with_control=False, binsize=50, unstranded=False, qval_cutoff=0.05, fold_change=[2.], min_clip_cov=0, pooling=False):
	"""DOCSTRING
	Args
	Returns
	"""
	# fetch the IP tag counts to gene regions
	if unique_only:
		interval_ip = \
			count_gene_read_tags(bam_dict['ubam.ip'], gene, is_unique=True, unstranded=unstranded)
	else:
		interval_ip = \
			count_gene_read_tags(bam_dict['ubam.ip'], gene, is_unique=True, unstranded=unstranded) + \
			count_gene_read_tags(bam_dict['mbam.ip'], gene, is_unique=False, unstranded=unstranded)
	
	# pooling reads from different replicate if provided multiple replicates
	# and poolinging is turned on
	if interval_ip.shape[0]>1 and pooling:
		interval_ip = np.apply_along_axis(np.sum, 0, interval_ip).reshape((1, interval_ip.shape[1]))
	
	# skip if there are not enough reads
	ip_sum = np.apply_along_axis(np.sum, 1, interval_ip)
	valid_ip_sample = np.where(ip_sum > min_clip_cov)[0]
	if len(valid_ip_sample)==0:
		#print "no reads"
		return ''
	interval_ip = interval_ip[valid_ip_sample,:]
		
	# fetch/construct the input tag counts
	if with_control:
		## count control tags if available
		if unique_only:
			interval_con = \
				count_gene_read_tags(bam_dict['ubam.con'], gene, is_unique=True, unstranded=unstranded)
		else:
			interval_con = \
				count_gene_read_tags(bam_dict['ubam.con'], gene, is_unique=True, unstranded=unstranded) + \
				count_gene_read_tags(bam_dict['mbam.con'], gene, is_unique=False, unstranded=unstranded)
	else:
		## otherwise, construct a uniform *fake* control
		interval_con = \
				np.ones((1, interval_ip.shape[1]))*np.sum(interval_ip)/interval_ip.shape[1]
		#interval_con = np.empty(interval_ip.shape)
		#for i in range(interval_ip.shape[0]):
		#	interval_con[i, ] = \
		#		np.ones((interval_ip.shape[1]))*np.sum(interval_ip[i,])/interval_ip.shape[1]
	
	# pooling reads from different replicate if provided multiple replicates
	# and poolinging is turned on
	if interval_con.shape[0]>1 and pooling:
		interval_con = np.apply_along_axis(np.sum, 0, interval_con).reshape((1, interval_con.shape[1]))

	
	# bin tag counts into bins
	intv_bin_ip = bin_interval_counts(interval_ip, binsize=binsize)
	intv_bin_con = bin_interval_counts(interval_con, binsize=binsize)
	
	# perform statistical test
	signal_val, binscore_adj, binscore = test_bin_negbinom(intv_bin_ip, intv_bin_con, with_control=with_control)
	
	# DO NOT USE
	#if with_control or intv_bin_ip.shape[0]>1:
	#	signal_val, binscore_adj = test_bin_negbinom(intv_bin_ip, intv_bin_con, with_control=with_control)
	#else:
	#	signal_val, binscore_adj = bin_test_alternatives.test_bin_fisher(intv_bin_ip, intv_bin_con)
	
	if signal_val is None:
		return ''
	
	# build human-readable outputs
	## "narrowPeak" format from 
	## https://genome.ucsc.edu/FAQ/FAQformat.html#format12
	## chr start end name 1000 strand signalValue pVal qVal peak
	narrowPeak_formatter = "%s\t%i\t%i\t%s\t1000\t%s\t%s\t%.3e\t%.3e\t.\n"
	BED = ''
	if len(fold_change)==1:
		lb = np.log(fold_change[0]) if with_control else fold_change[0]
		ub = np.inf
	else:
		assert fold_change[0]<fold_change[1]
		lb = np.log(fold_change[0]) if with_control else fold_change[0]
		ub = np.log(fold_change[1]) if with_control else fold_change[1]
	peak_num = 0
	for i in range(len(binscore_adj)):
		qval = binscore_adj[i]
		pval = binscore[i]
		signal = signal_val[i]
		num_signal = float(signal) if with_control else np.mean([float(x) for x in signal.split(',')])
		if qval < qval_cutoff and num_signal > lb and num_signal < ub :
			chr = gene[0]
			binstart = gene[1] + i*binsize
			binend = gene[1] + (i+1)*binsize
			strand = gene[3]
			peak_num += 1
			peak_name = gene[4] + '-%i'%peak_num
			BED += narrowPeak_formatter % (chr, binstart, binend, peak_name, strand, signal, pval, qval)
	return BED
	

def _child_peak_caller( (ip_bam_list, con_bam_list, child_gene_list, gene_annot, unique_only, with_control, unstranded, binsize, qval_cutoff, fold_change, min_clip_cov, pooling) ):
	"""DOCSTRING
	Args
	Returns
	"""
	# open file handler for child process
	bam_dict = make_bam_handler_dict(ip_bam_list, con_bam_list)
	
	# do the child's jobs
	BED = ''
	pid = os.getpid()
	tot = len(child_gene_list)
	for i in range(len(child_gene_list)):
		if not i%200:
			logger.debug('pid %s : %i / %i (%.2f%%)'% (pid, i, tot, float(i)/float(tot)*100))
		gene_name = child_gene_list[i]
		gene = gene_annot[gene_name]
		BED += call_gene_peak(bam_dict, gene, 
				unique_only=unique_only, with_control=with_control, 
				unstranded=unstranded, binsize=binsize, qval_cutoff=qval_cutoff, fold_change=fold_change,
				min_clip_cov=min_clip_cov, pooling=pooling)
	# close the handler
	_ = map(lambda x: x.close(), [bam for x in bam_dict.values() for bam in x])
	return BED


def make_bam_handler_dict(ip_bam_list, con_bam_list):
	"""DOCSTRING
	Args
	Returns
	"""
	bam_dict = defaultdict(list)
	try:
		ubam_ip, mbam_ip = ip_bam_list
	except:
		ubam_ip = ip_bam_list[0]
		mbam_ip = None
	for bam_fn in ubam_ip.split(','):
		if not os.path.isfile(bam_fn):
			raise Exception('%s not found'%bam_fn)
		bam_dict['ubam.ip'].append( pysam.Samfile(bam_fn, 'rb') )
	if mbam_ip is not None:
		for bam_fn in mbam_ip.split(','):
			if not os.path.isfile(bam_fn):
				raise Exception('%s not found'%bam_fn)
			bam_dict['mbam.ip'].append( pysam.Samfile(bam_fn, 'rb') )
		
	if con_bam_list is None:
		return bam_dict
	try:
		ubam_con, mbam_con = con_bam_list
	except:
		ubam_con = con_bam_list[0]
		mbam_con = None
	for bam_fn in ubam_con.split(','):
		if not os.path.isfile(bam_fn):
			raise Exception('%s not found'%bam_fn)
		bam_dict['ubam.con'].append( pysam.Samfile(bam_fn, 'rb') )
	if mbam_con is not None:
		for bam_fn in mbam_con.split(','):
			if not os.path.isfile(bam_fn):
				raise Exception('%s not found'%bam_fn)
			bam_dict['mbam.con'].append( pysam.Samfile(bam_fn, 'rb') )
	return bam_dict


def peakcaller(ip_bam_list, gtf_fp, con_bam_list=None, nthread=8, 
		out_dir='.', binsize=50,
		unique_only=False, unstranded=True,
		qval_cutoff=0.05, fold_change=[2.], 
		min_clip_cov=0, pooling=False):
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
	# indicator for whether we have real control
	if con_bam_list is not None:
		with_control = True
	else:
		with_control = False
	if len(ip_bam_list)==1 and not unique_only:
		logger.info('with only 1 bam file provided, i can only run --unique-only')
		unique_only=True

	if nthread == 1:
		# file handlers only for single-thread
		bam_dict = make_bam_handler_dict(ip_bam_list, con_bam_list)
	
	if unique_only:
		ofile = open(os.path.join(out_dir, 'narrow_peak.unique.bed'), 'w')
	else:
		ofile = open(os.path.join(out_dir, 'narrow_peak.combined.bed'), 'w')
	
	# read in GTF
	gene_annot = read_gtf(gtf_fp)
	
	if nthread == 1:
	# iteratively call peaks in each gene
		peak_counter = 0
		for gene_name in tqdm(gene_annot):
			gene = gene_annot[gene_name]
			BED = call_gene_peak(bam_dict, gene, 
				unique_only=unique_only, with_control=with_control, 
				unstranded=unstranded, binsize=binsize, qval_cutoff=qval_cutoff, fold_change=fold_change,
				min_clip_cov=min_clip_cov, pooling=pooling)
			ofile.write(BED)
			#print BED
			peak_counter += len(BED.split('\n'))-1
		
		# close the handler
		ofile.close()
		_ = map(lambda x: x.close(), [bam for x in bam_dict.values() for bam in x])	
	else:
	# multi-threading on subset of genes
		logger.info('multi-threading')
		gene_list = gene_annot.keys()
		child_gene_list = [x for x in chunkify(gene_list, nthread)]
		pool=Pool(processes=nthread)
		BED_list = pool.map(
			_child_peak_caller, 
			[(ip_bam_list, con_bam_list, child_gene_list[i], gene_annot, unique_only, with_control, unstranded, binsize, qval_cutoff, fold_change, min_clip_cov, pooling) for i in range(nthread)]
			)
		pool.terminate()
		pool.join()
		BED = ''.join(BED_list)
		peak_counter = len(BED.split('\n'))-1
		ofile.write(BED)
		ofile.close()
	logger.info('called %i peaks'%peak_counter)
	return
	

def chunkify(a, n):
	"""Separate a list (a) into consecutive n chunks.
	Args:
	Returns:
		the chunkified index
	"""
	k, m = len(a) / n, len(a) % n
	return (a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in xrange(n))


def parser(args):
	"""DOCSTRING
	Args
	Returns
	"""
	try:
		ip_bam_list = args.in_bam
		con_bam_list = args.con_bam
		out_dir = args.out_dir
		if not os.path.isdir(out_dir):
			os.mkdir(out_dir)
		gtf_fp = args.gtf_fp
		nthread = args.nthread
		unstranded = args.unstranded
		unique_only = args.unique_only
		binsize = args.binsize
		qval_cutoff = args.qval_cutoff
		fold_change = args.fold_change
		min_clip_cov = args.min_clip_cov
		pooling = args.pooling
		logger = logging.getLogger('CLAM.Peakcaller')
		logger.info('start')
		logger.info('run info: %s'%(' '.join(sys.argv)))
		
		peakcaller(ip_bam_list, gtf_fp, con_bam_list, nthread, 
			out_dir=out_dir, binsize=binsize,
			unique_only=unique_only, unstranded=unstranded,
			qval_cutoff=qval_cutoff,
			fold_change=fold_change,
			min_clip_cov=min_clip_cov,
			pooling=pooling)
		
		logger.info('end')
	except KeyboardInterrupt():
		sys.exit(0)
	return


if __name__ == '__main__':
	### set up logger
	logger = logging.getLogger('CLAM')
	logger.setLevel(logging.DEBUG)
	# create file handler which logs info messages
	fh = logging.FileHandler(
		'CLAM.Peakcaller.'+'-'.join(str(datetime.datetime.now()).replace(':','-').split()) + '.log')
	fh.setLevel(logging.INFO)
	# create console handler with even debug log level
	ch = logging.StreamHandler()
	ch.setLevel(logging.DEBUG)
	# create formatter and add it to the handlers
	formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s -\n %(message)s')
	fh.setFormatter(formatter)
	ch.setFormatter(formatter)
	# add the handlers to the logger
	logger.addHandler(fh)
	logger.addHandler(ch)
	###
	logger.info('start')
	logger.info('run info: %s'%(' '.join(sys.argv)))
	
	ip_bam_list, con_bam_list, unique_only, nthread = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
	if len(sys.argv)>5:
		out_dir = sys.argv[5]
	else:
		out_dir = '.'
	if len(sys.argv)>6:
		unstranded = True if sys.argv[5]!='0' else False
	else:
		unstranded = False
	unique_only = False if unique_only=='0' else True
	if con_bam_list == 'None':
		con_bam_list = None
	gtf_fp = '/u/nobackup/yxing/NOBACKUP/frankwoe/hg19/gencode.v19.annotation.gtf'
	peakcaller(ip_bam_list, gtf_fp, con_bam_list=con_bam_list, 
		unique_only=unique_only, nthread=nthread, out_dir=out_dir)
	logger.info('end')
