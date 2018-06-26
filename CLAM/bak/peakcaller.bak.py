#!/usr/bin/env python

"""
This peak-caller script is part of the CLAM pipeline.

It takes input from re-aligner output, and use permutation to call peaks.

Tested under python 2.7.3
"""

__author__ = 'Zijun Zhang'
__version__ = '1.1.0'
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
from stats import ztnb_em


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
				gene_id = re.search(r'gene_id "(.+?)"', ele[-1]).group(1)
			except AttributeError:
				continue
			gene_annot[gene_id] = [chr, start, end, strand]
	return gene_annot


def count_gene_read_tags(bam_list, (chr, start, end, strand), is_unique=True, unstranded=False):
	""" count the tagger positions for all reads in a given genomic interval
	Args:
	Returns:
	"""
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
			read_tags = [ (x.opt('RT'), x.opt('AS')) for x in bam.fetch(chr, start, end) \
				if unstranded or x.is_reverse==is_reverse]
		
		for tag in read_tags:
			if tag[0]<start or tag[0]>=end:
				continue
			interval[i, tag[0]-start] += tag[1]
	return interval


def	bin_interval_counts(interval, winsize=50):
	bins = np.zeros( ( interval.shape[0], int(np.ceil(interval.shape[1]/float(winsize))) ) )
	for i in range(bins.shape[1]):
		for j in range(interval.shape[0]):
			start, end = i*winsize, (i+1)*winsize-1
			bins[j, i] = np.sum(interval[j, start:end])
	return bins


def test_bin_negbinom(intv_bin_ip, intv_bin_con, correction_method='fdr_bh'):
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
	binsignal = np.empty(intv_counter)
	alpha_ip_vec = np.empty(intv_bin_ip.shape[0])
	alpha_con_vec = np.empty(intv_bin_con.shape[0])
	ip_sum = np.apply_along_axis(np.sum, 1, intv_bin_ip)
	con_sum = np.apply_along_axis(np.sum, 1, intv_bin_con)
	
	
	# compute the dispersion parameters
	for i in range(intv_bin_con.shape[0]):
		height = ztnb_em.collapse_data(intv_bin_con[i,])
		height[0] = 0
		ll, mu, alpha = ztnb_em.EM_estim_params(height, max_iter=100, verbose=False)
		alpha_con_vec[i] = alpha
	
	max_alpha = np.max(alpha_con_vec)
	for i in range(intv_bin_ip.shape[0]):
		height = ztnb_em.collapse_data(intv_bin_ip[i,])
		height[0] = 0
		ll, mu, alpha = ztnb_em.EM_estim_params(height, max_iter=100, verbose=False)
		alpha = max_alpha if alpha>max_alpha else alpha
		alpha_ip_vec[i] = alpha
	
	
	# perform test on each bin
	for i in range(intv_counter):
		this_ip = intv_bin_ip[:, i]
		others_ip = ip_sum - this_ip
		this_con = intv_bin_con[:, i]
		others_con = con_sum - this_con
		if np.sum(this_ip) == 0:
			binsignal[i], binscore[i] = np.nan, np.nan
			continue
		data = {
				'this_ip':np.round(this_ip),
				'others_ip':np.round(others_ip),
				'this_con':np.round(this_con),
				'others_con':np.round(others_con)
			}
		## constrained likelihood
		res_constrain = optimize.minimize(
				x0=np.ones(1+this_ip.shape[0]+others_ip.shape[0]), 
				fun=_neg_loglik_constrain,
				args=(data),
				method='bfgs',
				options={'disp':False}
			)
		## unconstrained likelihood
		res_unconstrain = optimize.minimize(
				x0=np.ones(2+this_ip.shape[0]+others_ip.shape[0]), 
				fun=_neg_loglik_unconstrain,
				args=(data),
				method='bfgs',
				options={'disp':False}
			)
		
		delta_mle = res_unconstrain.x[1]
		pval = 1 - chi2.cdf(2*(res_constrain.fun - res_unconstrain.fun), 1)
		binscore[i] = pval
		binsignal[i] = delta_mle
	
	# correcting for multiple-testing
	adj = multipletests(binscore[~ np.isnan(binscore)], alpha=0.05, method=correction_method)
	binscore_adj = np.asarray(binscore)
	binscore_adj[ ~ np.isnan(binscore) ] = adj[1]
	return binsignal, binscore_adj


def test_bin_poisson(intv_bin_ip, intv_bin_con, correction_method='fdr_bh'):
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
		
	def _neg_loglik_unconstrain(par, data):
		(l1, l2, l3, l4) = _par_to_vec(par, data, False)
		ll = np.sum(poisson.logpmf(data['this_ip'], mu=l1)) + \
			np.sum(poisson.logpmf(data['others_ip'], mu=l2)) + \
			np.sum(poisson.logpmf(data['this_con'], mu=l3)) + \
			np.sum(poisson.logpmf(data['others_con'], mu=l4))
		return -ll
	
	def _neg_loglik_constrain(par, data):
		(l1, l2, l3, l4) = _par_to_vec(par, data, True)
		ll = np.sum(poisson.logpmf(data['this_ip'], mu=l1)) + \
			np.sum(poisson.logpmf(data['others_ip'], mu=l2)) + \
			np.sum(poisson.logpmf(data['this_con'], mu=l3)) + \
			np.sum(poisson.logpmf(data['others_con'], mu=l4))
		return -ll
		
	intv_counter = intv_bin_ip.shape[1]
	assert intv_counter == intv_bin_con.shape[1]
	binscore = np.empty(intv_counter)
	binsignal = np.empty(intv_counter)
	ip_sum = np.apply_along_axis(np.sum, 1, intv_bin_ip)
	con_sum = np.apply_along_axis(np.sum, 1, intv_bin_con)
	for i in range(intv_counter):
		this_ip = intv_bin_ip[:, i]
		others_ip = ip_sum - this_ip
		this_con = intv_bin_con[:, i]
		others_con = con_sum - this_con
		if this_ip == 0:
			binsignal[i], binscore[i] = np.nan, 1.0
			continue
		## because Poisson (and other count-based methods) only
		## takes integers, here we take the floor of the fractional
		## multi-reads as a conservative approach
		data = {
				'this_ip':np.floor(this_ip),
				'others_ip':np.floor(others_ip),
				'this_con':np.floor(this_con),
				'others_con':np.floor(others_con)
			}
		
		res_constrain = optimize.minimize(
				x0=np.ones(1+this_ip.shape[0]+others_ip.shape[0]), 
				fun=_neg_loglik_constrain,
				args=(data),
				method='Nelder-Mead',
				options={'disp':False}
			)
		
		res_unconstrain = optimize.minimize(
				x0=np.ones(2+this_ip.shape[0]+others_ip.shape[0]), 
				fun=_neg_loglik_unconstrain,
				args=(data),
				method='Nelder-Mead',
				options={'disp':False}
			)
		
		delta_mle = res_unconstrain.x[1]
		pval = 1 - chi2.cdf(2*(res_constrain.fun - res_unconstrain.fun), 1)
		binscore[i] = pval
		binsignal[i] = delta_mle
	adj = multipletests(binscore, alpha=0.05, method=correction_method)
	binscore_adj = adj[1]
	return binsignal, binscore_adj


def test_bin_fisher(intv_bin_ip, intv_bin_con, with_control=True, correction_method='fdr_bh'):
	"""DOCSTRING
	Args
	Returns
	"""
	if intv_bin_ip.shape[0] != 1:
		raise Exception('Fisher exact test does not deal with replicates.')
	intv_counter = intv_bin_ip.shape[1]
	assert intv_counter == intv_bin_con.shape[1]
	binscore = np.empty(intv_counter)
	binsignal = np.empty(intv_counter)
	ip_sum = np.sum(intv_bin_ip[0,])
	con_sum = np.sum(intv_bin_con[0,])
	for i in range(intv_counter):
		this_ip = intv_bin_ip[0, i]
		others_ip = ip_sum - this_ip
		this_con = intv_bin_con[0, i]
		others_con = con_sum - this_con
		if this_ip == 0:
			binsignal[i], binscore[i] = np.nan, 1.0
			continue
		_, binscore[i] = fisher_exact([[this_ip, others_ip], [this_con, others_con]], alternative='greater')
		if with_control:
			binsignal[i] = this_ip/others_ip / this_con*others_con
		else:
			binsignal[i] = this_ip
		
	adj = multipletests(binscore, alpha=0.05, method=correction_method)
	binscore_adj = adj[1]
	return binsignal, binscore_adj


def call_gene_peak(bam_dict, gene, unique_only=False, with_control=False, winsize=50, unstranded=False):
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
	
	# skip if there are no reads
	if np.sum(interval_ip) == 0:
		#print "no reads"
		return ''
		
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
	
	# bin tag counts into bins
	intv_bin_ip = bin_interval_counts(interval_ip, winsize=winsize)
	intv_bin_con = bin_interval_counts(interval_con, winsize=winsize)
	
	# perform statistical test
	signal_val, binscore_adj = test_bin_negbinom(intv_bin_ip, intv_bin_con)
	#signal_val, binscore_adj = test_bin_poisson(intv_bin_ip, intv_bin_con)
	#signal_val, binscore_adj = test_bin_fisher(bin_interval_ip, bin_interval_con, with_control=with_control)
	
	# build human-readable outputs
	## "narrowPeak" format from 
	## https://genome.ucsc.edu/FAQ/FAQformat.html#format12
	## chr start end name 1000 strand signalValue pVal qVal peak
	narrowPeak_formatter = "%s\t%i\t%i\t.\t1000\t%s\t%s\t.\t%.3f\t.\n"
	BED = ''
	for i in range(len(binscore_adj)):
		qval = binscore_adj[i]
		signal = signal_val[i]
		if qval<0.05:
			chr = gene[0]
			binstart = gene[1] + i*winsize
			binend = gene[1] + (i+1)*winsize-1
			strand = gene[3]
			BED += narrowPeak_formatter % (chr, binstart, binend, strand, signal, qval)
	return BED
	


def peakcaller(tmp_dir, out_dir, gtf_fp, unique_only=False, with_replicates=False, with_control=False, unstranded=False):
	"""DOCSTRING
	Args:
	Returns:
	"""
	# file handlers
	mbam = pysam.Samfile(os.path.join(out_dir, 'realigned.sorted.bam'),'rb')
	ubam = pysam.Samfile(os.path.join(tmp_dir, 'unique.sorted.bam'),'rb')
	bam_dict = {'ubam.ip':[ubam], 'mbam.ip':[mbam]}
	if unique_only:
		ofile = open(os.path.join(out_dir, 'narrow_peaks.unique.bed'), 'w')
	else:
		ofile = open(os.path.join(out_dir, 'narrow_peaks.combined.bed'), 'w')
	
	# read in GTF
	gene_annot = read_gtf(gtf_fp)
	
	# iteratively call peaks in each gene
	peak_counter = 0
	for gene_name in tqdm(gene_annot):
		gene = gene_annot[gene_name]
		BED = call_gene_peak(bam_dict, gene, 
			unique_only=unique_only, with_control=with_control, 
			unstranded=unstranded)
		ofile.write(BED)
		#print BED
		peak_counter += len(BED.split('\n'))
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


if __name__ == '__main__':
	### set up logger
	logger = logging.getLogger('CLAM')
	logger.setLevel(logging.DEBUG)
	# create file handler which logs even debug messages
	fh = logging.FileHandler(
		'CLAM.Peakcaller.'+'-'.join(str(datetime.datetime.now()).replace(':','-').split()) + '.log')
	fh.setLevel(logging.DEBUG)
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
	###
	logger.info('start')
	logger.info('run info: %s'%(' '.join(sys.argv)))
	
	tmp_dir, out_dir, unique_only = sys.argv[1], sys.argv[2], sys.argv[3]
	unique_only = False if unique_only=='0' else True
	gtf_fp = '/u/nobackup/yxing/NOBACKUP/frankwoe/hg19/gencode.v19.annotation.gtf'
	peakcaller(tmp_dir, out_dir, gtf_fp, unique_only=unique_only)
	logger.info('end')
