#!/usr/bin/python

"""
This peak-caller script is part of the CLAM pipeline.

It takes input from re-aligner output, and use permutation to call peaks.

Tested under python 2.7.3
"""

__author__ = 'Zijun Zhang'
__version__ = '1.0.0'
__email__ = 'zj.z@ucla.edu'


from optparse import OptionParser
import os, subprocess, sys
from collections import defaultdict
from statsmodels.sandbox.stats.multicomp import multipletests
from time import strftime
import cPickle as pickle
import bisect, random
import pysam
import pybedtools
from multiprocessing import Pool

def main():
	"""
	The main wrapper for CLAM peak-caller.
	"""
	# options parsing
	usage='usage: %prog <options>'
	parser=OptionParser(usage)
	
	parser.add_option('--resume', dest='resume', action='store_true', default=False, help='Resume mode - skipping pre-processing [Default: %default]')
	parser.add_option('--verbose', dest='verbose', action='store_true', default=False, help='Verbose mode - print out all intermediate steps [Default: %default]')
	parser.add_option('-o', dest='output_dir', default='./out_CLAM', help='Output file folder [Default %default]')
	parser.add_option('-t', dest='tmp_dir', default='./tmp_CLAM', help='Temporary file folder [Default %default]')
	parser.add_option('-p', dest='peak_file', default=None, help='Output peak calling filename; if None then do not call peaks  [Default %default]')
	parser.add_option('--is-stranded', dest='is_stranded', default=False, action='store_true', help='Indicates if the reads are mapped with strand information. [Default: %default]')
	parser.add_option('--extend', dest='extend', type='int', default=50, help='Extend to given nucleotides symmetrically at peak calling [Default: %default]')
	parser.add_option('--pval-cutoff', dest='pval_cutoff', type='float', default=0.001, help='Corrected p-value threshold at peak calling [Default: %default]')
	parser.add_option('--merge-size', dest='merge_size', type='int', default=50, help='merging window size at peak calling [Default: %default]')
	parser.add_option('--max-iter', dest='max_iter', type='int', default=1000, help='maximum iterations for permutation tests [Default: %default]')
	parser.add_option('-g', dest='gtf', default='./GTF/hg19_ensembl.sorted_gene.bed', help='GTF file [Default: %default]')
	parser.add_option('--ThreadN', dest='nb_proc', type='int', default=4, help='Number of threads when doing permutations. [Default: %default]')
	parser.add_option('--seed', dest='seed', type='int', default=100, help='Random seed for permutations. [Default: %default]')
	parser.add_option('--merge-method', dest='merge_method', type='int', default=1, help='Peak merging method. 1: Single Nucleotide 2: Broad Coverage [Default: %default]')
	parser.add_option('--pval-method', dest='correction_method', type='int', default=1, help='Multiple testing correction method. 1: Bonferroni 2: BH FDR [Default: %default]')
	parser.add_option('--call-transcriptome', dest='call_all', action='store_true', default=False, help='Call peaks on transcriptome instead of genes with multi-mappers. [Default: %default]')
	
	(options,args)=parser.parse_args()
	
	output_dir=os.path.abspath(options.output_dir)
	tmp_dir=os.path.abspath(options.tmp_dir)
	verbose=options.verbose
	
	#random.seed(options.seed)
	
	write_parameter_log(options, output_dir)
	
	# find transcripts associated with multi-mapped reads
	if verbose:
		print_time_stamp('Finding transcripts with multimapped reads.')
	if not os.path.isfile(output_dir + '/CLAM_mapper.sorted.out'):
		subprocess.call(''' sort -k1,1 -k2,2n %s/CLAM_mapper.out | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' > %s/CLAM_mapper.sorted.out ''' % (output_dir, output_dir), shell=True)
	# Note: tid_list: tid -> [chr:strand, start, end]
	tid_list=read_aligner_output(output_dir + '/CLAM_mapper.sorted.out', options.gtf, options.is_stranded, tmp_dir, options.resume, options.call_all)
	
	# make bam file for re-aligner output, if non-exist
	if not (options.resume and os.path.isfile(output_dir + '/assigned_multimapped_reads.bam')):
		if verbose:
			print_time_stamp('Making bamfile for aligner output.')
		header_cmd='samtools view -H ' + tmp_dir + '/filter100.sorted.bam > ' + output_dir + '/sam_header.sam'
		subprocess.call(header_cmd, shell=True)
		body_cmd = ''' awk '{if($6=="+"){print $4"\t256\t"$1"\t"$2+1"\t0\t"$3-$2+1"M\t*\t0\t0\t*\t*\tAS:f:"$5}else{print $4"\t272\t"$1"\t"$2+1"\t0\t"$3-$2+1"M\t*\t0\t0\t*\t*\tAS:f:"$5 }}' ''' + output_dir + '/CLAM_mapper.sorted.out > ' + output_dir + '/CLAM_mapper.sorted.sam'
		subprocess.call(body_cmd, shell=True)
		makeBam_cmd = 'cat %s/sam_header.sam %s/CLAM_mapper.sorted.sam | samtools view -bS - > %s/assigned_multimapped_reads.bam' % (output_dir, output_dir,output_dir)
		subprocess.call(makeBam_cmd, shell=True)
		index_cmd = 'samtools index %s/assigned_multimapped_reads.bam' % output_dir
		subprocess.call(index_cmd, shell=True)
	
	# multi-processing peak-caller
	if not (options.resume and os.path.isfile(tmp_dir+'/unique_to_qval.pdata') and os.path.isfile(tmp_dir+'/combined_to_qval.pdata')):
		child_transcr_ind = list(chunkify(range(len(tid_list)), options.nb_proc))
		
		pool = Pool(processes=options.nb_proc)
		
		unibam_file=tmp_dir+'/filter100.sorted.bam'
		multibam_file=output_dir+'/assigned_multimapped_reads.bam'
		tid_to_qval_compact = pool.map(get_permutation_fdr, [ (unibam_file, multibam_file, tid_list, child_transcr_ind[i], options.pval_cutoff, options.max_iter, options.is_stranded, verbose, options.correction_method, options.seed) for i in range(options.nb_proc) ])
		
		pool.terminate()
		pool.join()
		
		unique_tid_to_qval, combined_tid_to_qval = unpack_tid_to_qval(tid_to_qval_compact)
		pickle.dump(unique_tid_to_qval, open(tmp_dir+'/unique_to_qval.pdata','wb'), -1)
		pickle.dump(combined_tid_to_qval, open(tmp_dir+'/combined_to_qval.pdata','wb'), -1)
	else:
		print_time_stamp('Resume mode, found qval data files.')
		unique_tid_to_qval=pickle.load(open(tmp_dir+'/unique_to_qval.pdata','rb'))
		combined_tid_to_qval=pickle.load(open(tmp_dir+'/combined_to_qval.pdata','rb'))
	
	# merge peaks
	if options.merge_method==1:
		merge_peaks=merge_peaks_singleNucl
		mm='singleNucl'
	elif options.merge_method==2:
		merge_peaks=merge_peaks_broadPeak
		mm='broadPeak'
	else:
		merge_peaks=merge_peaks_singleNucl
		mm='unknown selection, using default singleNucl'
		
	if verbose:
		print_time_stamp('Merging peaks within ' + str(options.merge_size) + 'bp, using ' + mm + '..')
	
	unique_peaks=merge_peaks(unique_tid_to_qval, options.merge_size, options.pval_cutoff)
	combined_peaks=merge_peaks(combined_tid_to_qval, options.merge_size, options.pval_cutoff)
	
	print_time_stamp('Comparing results and writing to file..')
	
	# write peak-calling results to file.
	with open(output_dir + '/all_peaks.txt', 'w') as f:
		for peak in unique_peaks:  # peak = ['chr\tstart\tend\tstrand', 'height\tqval\t', tid]
			if options.extend is None:
				wt_loc=peak[0]
			else:
				wt_loc=extend_peak_region(peak[0], options.extend)
			f.write(wt_loc + '\t' + '\t'.join([str(x) for x in peak[1]]) + '\t' + peak[2] + '\tunique\n')
		for peak in combined_peaks:
			if options.extend is None:
				wt_loc=peak[0]
			else:
				wt_loc=extend_peak_region(peak[0], options.extend)
			f.write(wt_loc + '\t' + '\t'.join([str(x) for x in peak[1]]) + '\t' + peak[2] + '\tcombined\n')
	subprocess.call(''' sort -k1,1 -k2,2n %s/all_peaks.txt | awk '{print $1"\t"$2"\t"$3"\t"$5";"$6";"$7"\t"$8"\t"$4}'  | bedtools merge -s -d -1 -i stdin -c 4,5,6, -o collapse,collapse,distinct  > %s''' % (output_dir, options.peak_file), shell=True)
		
	print_time_stamp('Peak-calling done.')

def write_parameter_log(options, output_dir):
	"""
	Write paramter values to a log file, named by current time.
	"""
	merge_method_dict={1:'singleNucl', 2:'broadPeak'}
	correction_method_dict={1:'Bonferroni', 2:'BH_FDR'}
	with open(output_dir+'/CLAM_Peaker.Parameters.'+ strftime("%Y%m%d_%H%M") + '.txt', 'w') as log:
		log.write('CLAM Peaker ' + __version__ + '\n')
		log.write('resume: ' + str(options.resume) + '\n')
		log.write('verbose: ' + str(options.verbose) + '\n')
		log.write('output_dir:' + str(options.output_dir) + '\n')
		log.write('tmp_dir: ' + str(options.tmp_dir) + '\n')
		log.write('peak_file: ' + str(options.peak_file) + '\n')
		log.write('is_stranded: ' + str(options.is_stranded) + '\n')
		log.write('extend: ' + str(options.extend) + '\n')
		log.write('pval_cutoff: ' + str(options.pval_cutoff) + '\n')
		log.write('merge_size: ' + str(options.merge_size) + '\n')
		log.write('max_iter: ' + str(options.max_iter) + '\n')
		log.write('gtf: ' + str(options.gtf) + '\n')
		log.write('seed: ' + str(options.seed) + '\n')
		log.write('merge_method: ' + merge_method_dict[options.merge_method] + '\n')
		log.write('correction_method: ' + correction_method_dict[options.correction_method] + '\n')
		log.write('thread: ' + str(options.nb_proc) + '\n')

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
	
def get_permutation_fdr((unibam_file, multibam_file, tid_list, tid_ind, pval_cutoff, max_iter, is_stranded, verbose, correction_method, seed)):
	"""
	General permutation wrapper for a list of genes. Gets called by multi-processing generated by Pool()
	Returns packed FDRs from each child process.
	"""
	random.seed(seed)
	
	unique_tid_to_qval=defaultdict(list)
	combined_tid_to_qval=defaultdict(list)
	
	unibam=pysam.Samfile(unibam_file, 'rb')
	multibam=pysam.Samfile(multibam_file, 'rb')
	
	processed=0
	pid=os.getpid()
	
	for ind in tid_ind:
		processed+=1
		if verbose and not processed % 100:
			print_time_stamp(str(processed) + '/' + str(len(tid_ind)) + ' finished in pid ' + str(pid))
		tid, chr, strand, start, end = tid_list[ind]
		unique_reads = read_tid_frag_from_bam(tid_list[ind], unibam, is_stranded, True)
		multi_reads = read_tid_frag_from_bam(tid_list[ind], multibam, is_stranded, False)
		
		this_unique_to_qval = do_permutation(tid_list[ind], unique_reads, max_iter, pval_cutoff, correction_method)
		this_combined_to_qval = do_permutation(tid_list[ind], unique_reads+multi_reads, max_iter, pval_cutoff, correction_method)
		
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
	tid, chr, strand, tstart, tend = transcr
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

	if correction_method==2:
		qval_list=multipletests(pval_list, method='fdr_bh')[1]
	else:
		qval_list=[min(x*(len(set([int(y) for y in height_to_pval if y!=0]))), 1.0) for x in pval_list]
	
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

def read_aligner_output(rm_out, gtffile, is_stranded, tmp_dir, resume, call_all):
	"""
	Use bedtools to get transcripts/genes with multi-mapped reads.
	Returns a list of transcripts/genes.
	"""
	if not (resume and os.path.isfile(tmp_dir + '/gtf2multireads.bed')):
		rm_bed=pybedtools.BedTool(rm_out)
		gtf=pybedtools.BedTool(gtffile)
		gtf_bed_rm = gtf.intersect(rm_bed, s=True, u=True) if is_stranded else gtf.intersect(rm_bed, u=True)
		gtf_bed_rm.saveas(tmp_dir + '/gtf2multireads.bed')
		pybedtools.cleanup()
	
	tid_list=[]
	if call_all:
		gtf_to_read=gtffile
	else:
		gtf_to_read=tmp_dir+'/gtf2multireads.bed'
	with open(gtf_to_read,'r') as f:
		for line in f:
			ele=line.rstrip().split('\t')
			gene_id=ele[3]
			gene_chr, gene_start, gene_end=ele[0], int(ele[1]), int(ele[2])
			gene_strand=ele[5]
			tid_list.append([gene_id, gene_chr, gene_strand, gene_start, gene_end])
	print_time_stamp('Read transcripts with multi-reads: ' + str(len(tid_list)))
	return tid_list

def read_tid_frag_from_bam(tid, bamfile, is_stranded, is_unique):
	"""
	Use pysam to fetch reads info for a given gene and its loci.
	Returns reads, read weights and its mapped loci.
	"""
	tid_reads=[]
	gene, chr, strand, start, end=tid
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
		read_length = read.qlen if read.qlen > 0 else read.positions[-1] - read.positions[0] + 1
		if read.pos-start>=0 and read_length<500:    # to avoid junction reads
			tid_reads.append([read.qname, read.pos-start, read_length, score])
	return tid_reads

def print_time_stamp(msg):
	"""
	Reporter function for logging.
	"""
	current_time='[' + strftime("%Y-%m-%d %H:%M:%S") + '] '
	print >> sys.stderr, current_time + msg

if __name__=='__main__':
	main()