""" parseBAM from Yan Gao
https://github.com/yangao07/pyParseBAM/blob/master/parse_bam.py
Author: Yan Gao
Date: 1.9.2018
"""

import sys, re
import pysam as ps
import utils as ut

#### cigar operation:
BAM_CMATCH = 0  # M
BAM_CINS = 1  # I
BAM_CDEL = 2  # D
BAM_CREF_SKIP = 3  # N
BAM_CSOFT_CLIP = 4  # S
BAM_CHARD_CLIP = 5  # H
BAM_CPAD = 6  # P
BAM_CEQUAL = 7  # =
BAM_CDIFF = 8  # X
BAM_CBACK = 9  # B


# cigar stats:
# M	BAM_CMATCH	0
# I	BAM_CINS	1
# D	BAM_CDEL	2
# N	BAM_CREF_SKIP	3
# S	BAM_CSOFT_CLIP	4
# H	BAM_CHARD_CLIP	5
# P	BAM_CPAD	6
# =	BAM_CEQUAL	7
# X	BAM_CDIFF	8
# B	BAM_CBACK	9
# NM	NM tag	10
def get_ref_op_length(cigar_stats=[]):
    # get op length for MDNP=X
    op_len = cigar_stats[0] + cigar_stats[2] + cigar_stats[3] + cigar_stats[6] + cigar_stats[7] + cigar_stats[8]
    return op_len


def get_read_op_length(cigar_stats):
    # get op length for MISH=X
    op_len = cigar_stats[0] + cigar_stats[1] + cigar_stats[4] + cigar_stats[5] + cigar_stats[7] + cigar_stats[8]
    return op_len


def get_aligned_read_length(cigar_stats):
    # get op length for MI=X
    op_len = cigar_stats[0] + cigar_stats[1] + cigar_stats[7] + cigar_stats[8]
    return op_len


def minipulate_cigar(r=ps.AlignedSegment, old='', new=''):
    r.cigarstring = re.sub(r'%s' % old, new, r.cigarstring)


def get_spec_MD(mdstr='', start=0, end=0):  # '23AC20T', 10, 30 => ' 13AC5
    mSub = re.sub(r'([ACGTNacgtn])', ' \\1 ', mdstr)
    mSplit = re.split('[ ]+', mSub)
    start_remain_len = start
    end_remain_len = end - start
    ret_md = []
    mi = 0
    # print mSplit, start_remain_len, end_remain_len
    while start_remain_len > 0:
        if mSplit[mi].isdigit():
            if int(mSplit[mi]) > start_remain_len:
                mSplit[mi] = str(int(mSplit[mi]) - start_remain_len)
                start_remain_len = 0
                break
            else:
                start_remain_len -= int(mSplit[mi])
        else:  # isalpha()
            start_remain_len -= 1
        mi += 1
    while end_remain_len > 0:
        if mSplit[mi].isdigit():
            if int(mSplit[mi]) >= end_remain_len:
                ret_md.append(str(end_remain_len))
                end_remain_len = 0
                break
            else:
                end_remain_len -= int(mSplit[mi])
                ret_md.append(mSplit[mi])
        else:  # isalpha()
            end_remain_len -= 1
            ret_md.append(mSplit[mi])
        mi += 1
    # print ret_md
    return ret_md


# MISMATCH: read_pos(first), ref_pos(first), len, read_base, ref_base
# INSERTION: ins_read_pos(first), ins_ref_pos(left),  len, ins_base
# DELETION: del_read_pos(left), del_ref_pos(first), len, del_base
def get_error_from_MD(cigartuples=[], mdstr='', full_query_seq='', ref_start=0):
    mis, ins, dele = [], [], []
    last_error = ''
    md_i, m_pos = 0, 0
    mdSub = re.sub(r'([\\^][ACGTNacgtn]+)[0]*', ' \\1 ', mdstr)
    mdSplit = mdSub.rsplit()
    ref_pos, query_pos = ref_start, 0

    for tuples in cigartuples:
        if tuples[0] == BAM_CMATCH:
            m = mdSplit[md_i]

            if m.startswith('^'):
                ut.format_time(sys.stderr, 'get_error_from_MD', 'Unexpected MD string: {}\n'.format(mdstr))
                sys.exit(1)
            mSub = re.sub(r'([ACGTNacgtn])', ' \\1 ', m)
            m_len = sum(map(int, (re.sub(r'([ACGTNacgtn])', '1', mSub)).rsplit()))

            # from m_pos to m_pos + tuples[1]
            sub_ms = get_spec_MD(m, m_pos, m_pos + tuples[1])
            for ms in sub_ms:
                if ms.isalpha():  # MISMATCH
                    if last_error != 'MIS' or mis[-1][0] != query_pos - 1:
                        mis_error = [query_pos, ref_pos, 1, full_query_seq[query_pos], ms]
                        mis.append(mis_error)
                    else:  # last_error == 'MIS' and  mis[-1][1] == ap[0] - 1:
                        mis[-1][-3] += 1
                        mis[-1][-2] += full_query_seq[query_pos]
                        mis[-1][-1] += ms
                    query_pos += 1
                    ref_pos += 1
                    last_error = 'MIS'
                elif ms.isdigit():  # MATCH
                    query_pos += int(ms)
                    ref_pos += int(ms)

            if m_pos + tuples[1] == m_len:
                md_i += 1
                m_pos = 0
            elif m_pos + tuples[1] < m_len:
                m_pos += tuples[1]
            else:  #
                ut.format_time(sys.stderr, 'get_error_from_MD', 'Unexpected MD string: {}\n'.format(mdstr))
                sys.exit(1)
        elif tuples[0] == BAM_CDEL:
            m = mdSplit[md_i]
            if not m.startswith('^'):
                ut.format_time(sys.stderr, 'get_error_from_MD', 'Unexpected MD string: {}\n'.format(mdstr))
                sys.exit(1)
            del_error = [query_pos - 1, ref_pos, tuples[1], m[1:]]
            dele.append(del_error)
            ref_pos += tuples[1]
            last_error = 'DEL'
            md_i += 1
        elif tuples[0] == BAM_CINS:
            ins_error = [query_pos, ref_pos - 1, tuples[1], full_query_seq[query_pos:query_pos + tuples[1]]]
            ins.append(ins_error)
            query_pos += tuples[1]
            last_error = 'INS'
        elif tuples[0] == BAM_CSOFT_CLIP or tuples[0] == BAM_CHARD_CLIP:
            query_pos += tuples[1]
        elif tuples[0] == BAM_CREF_SKIP:
            ref_pos += tuples[1]
        else:
            ut.format_time(sys.stderr, 'get_error_from_MD', 'Unexpected cigar: {}\n'.format(cigartuples))
            sys.exit(1)

    return mis, ins, dele


def get_error_from_Cigar(cigartuples=[], full_query_seq='', align_ref_seq='', ref_start=0):
    mis, ins, dele = [], [], []
    last_error = ''
    ref_pos, query_pos = ref_start, 0
    for tuples in cigartuples:
        if tuples[0] == BAM_CMATCH:
            for q, r in zip(full_query_seq[query_pos:query_pos + tuples[1]],
                            align_ref_seq[ref_pos - ref_start:ref_pos - ref_start + tuples[1]]):
                if q != r:  # MISMATCH
                    if last_error != 'MIS' or mis[-1][0] != query_pos - 1:
                        mis_error = [query_pos, ref_pos, 1, q, r]
                        mis.append(mis_error)
                    else:  # last_error == 'MIS' and  mis[-1][1] == ap[0] - 1:
                        mis[-1][-3] += 1
                        mis[-1][-2] += q
                        mis[-1][-1] += r
                    last_error = 'MIS'
                ref_pos += 1
                query_pos += 1
        elif tuples[0] == BAM_CDEL:
            del_error = [query_pos - 1, ref_pos, tuples[1],
                         align_ref_seq[ref_pos - ref_start:ref_pos - ref_start + tuples[1]]]
            dele.append(del_error)
            last_error = 'DEL'
            ref_pos += tuples[1]
        elif tuples[0] == BAM_CINS:
            ins_error = [query_pos, ref_pos - 1, tuples[1], full_query_seq[query_pos:query_pos + tuples[1]]]
            ins.append(ins_error)
            last_error = 'INS'
            query_pos += tuples[1]
        elif tuples[0] == BAM_CHARD_CLIP or tuples[0] == BAM_CSOFT_CLIP:
            query_pos += tuples[1]
        elif tuples[0] == BAM_CREF_SKIP:
            ref_pos += tuples[1]
        else:
            ut.format_time(sys.stderr, 'get_error_from_Cigar', 'Unexpected cigar: {}\n'.format(cigartuples))
            sys.exit(1)

    return mis, ins, dele


def get_align_block(cigartuples=[], ref_start=0):
    align_block = []
    start = ref_start
    end = ref_start - 1
    for tuples in cigartuples:
        if tuples[0] == BAM_CMATCH or tuples[0] == BAM_CDEL:
            end += tuples[1]
        elif tuples[0] == BAM_CREF_SKIP:
            align_block.append((start, end))
            start = end + tuples[1] + 1
            end = start -1
    align_block.append((start, end))

    return align_block
	

# MISMATCH: read_pos(first), ref_pos(first), len, read_base, ref_base
# INSERTION: ins_read_pos(first), ins_ref_pos(left),  len, ins_base
# DELETION: del_read_pos(left), del_ref_pos(first), len, del_base
def get_align_detial(r, ref_fa):
    if r.cigartuples[0][0] == BAM_CSOFT_CLIP or r.cigartuples[0][0] == BAM_CHARD_CLIP:
        left_clip = r.cigartuples[0][1]
    else:
        left_clip = 0
    if r.cigartuples[-1][0] == BAM_CSOFT_CLIP or r.cigartuples[-1][0] == BAM_CHARD_CLIP:
        right_clip = r.cigartuples[-1][1]
    else:
        right_clip = 0

    if r.has_tag('MD'):
        mdstr = r.get_tag('MD')
        mis_err, ins_err, dele_err = get_error_from_MD(r.cigartuples, mdstr, r.query_sequence, r.reference_start)
    else:
        ref = ps.FastaFile(ref_fa)
        align_ref_seq = ref.fetch(r.reference_name, r.reference_start, r.reference_start + r.reference_length)
        mis_err, ins_err, dele_err = get_error_from_Cigar(r.cigartuples, r.query_sequence, align_ref_seq,
                                                             r.reference_start)
    return [r.is_reverse, r.infer_read_length(), r.reference_start, r.reference_length, left_clip, right_clip, mis_err,
            ins_err, dele_err]