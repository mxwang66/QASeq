#!/usr/bin/env python3
# -*- coding: utf-8 -*-


__author__ = 'Michael X. Wang'


"""
QASeq analysis code for short/long amplicon panels. 
"""


import os
import tarfile
import gzip
import pickle
import json
import pysam
import math
import multiprocessing
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.backends.backend_pdf # for saving figures to PDF
from matplotlib import pyplot as plt
from tqdm.notebook import tqdm
from scipy.optimize import curve_fit
from Bio import SeqIO
#from Bio import pairwise2 # deprecated
from Bio.Align import PairwiseAligner
from time import time
from concurrent import futures # create child process to flush RAM
from collections import Counter


# constants for adapter trim
MIN_AMP_LEN = 60
SEARCH_LEN = 10
SPACER_LEN = 4
UMI_LEN = 15

UMI_SPACER_LEN = UMI_LEN + SPACER_LEN
SPACER_SEARCH_LEN = SPACER_LEN + SEARCH_LEN
UMI_SEARCH_LEN = UMI_LEN + SEARCH_LEN
UMI_SPACER_SEARCH_LEN = UMI_SPACER_LEN + SEARCH_LEN

# for breast cancer 179-plex and 226-plex
MIN_AMP_LEN_LEGACY = 45

# constants for UMI analysis
PRIMER_VAL_LEN = 5 # number of primer validation bases
STATIC_TH = 3 # fixed family size threshold
TOP_N = 3 # top N largest UMI for dynamic cutoff
DYNAMIC_RATE = 0.05
UMI_COUNT_TH = 6

# constants for dowmsampling
DS_RATE_LS = [0.95, 0.9, 0.85, 0.8, 0.75, 0.7] # list of downsampling rate
CV_MIN_UMI_TH = 20 # minimum UMI count to calculate CV


def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))


def make_tarfile(source_dir, compress_type='gz', compresslevel=6, keep_source=True):
    """
    compress a directory
    :param source_dir: the directory to compress
    :param compress_type: 'gz', 'xz' or 'bz2'
    :param compresslevel: int, 1-9
    :param keep_source: True to keep the source directory
    """
    output_filename = source_dir + '.tar.' + compress_type
    with tarfile.open(output_filename, 'w:%s' % compress_type, compresslevel=compresslevel) as tar:
        tar.add(source_dir, arcname=os.path.basename(source_dir))

    if not keep_source:
        # delete source directory
        os.system('rm -r %s' % source_dir)


def revcomp(seq: str):
    """
    Reverse complement of DNA sequence
    """
    complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N', 
                  'a':'t', 't':'a', 'c':'g', 'g':'c', 'n':'n'}
    bases = [complement[b] for b in seq]
    bases.reverse()
    return ''.join(bases)


def seq2complex(seq_str: str):
    """
    Convert string sequence to a list of complex
    """
    dic = {'A':1, 'T':-1, 'C':1j, 'G':-1j, 'N': 0, 'a':1, 't':-1, 'c':1j, 'g':-1j, 'n': 0}
    seq_list = [dic[s] for s in list(seq_str)]
    return seq_list


def complex2seq(seq_complex: list):
    """
    Convert a list of complex back to string
    """
    seq_str = []
    for c in seq_complex:
        if c.real > c.imag:
            if -c.real < c.imag:
                seq_str.append('A')
            elif -c.real > c.imag:
                seq_str.append('G')
            else:
                seq_str.append('N')
        elif c.real < c.imag:
            if -c.real < c.imag:
                seq_str.append('C')
            elif -c.real > c.imag:
                seq_str.append('T')
            else:
                seq_str.append('N')
        else:
            seq_str.append('N')
    return ''.join(seq_str)


def variant_call(ref: str, var: str, n_chr='?', chr_pos=1):
    """
    find variants between a reference sequence and a instersted sequence
    HGVS Nomenclature is used in outputs
    :param ref: reference sequence
    :param var: sequence of interest
    :param n_chr: chromosome number of reference
    :param chr_pos: start of coordinate of reference on chromosome
    :return: a DataFrame of all the variants
    """
    # set gap penalty higher than match score, gap continue penalty equal to match score, and no mismatch penalty
    #align = pairwise2.align.globalms(ref, var, 2, 0, -3, -2, one_alignment_only=True)[0]
    aligner = PairwiseAligner(mode='global', match_score=2, mismatch_score=0, open_gap_score=-3, extend_gap_score=-2)
    alignment = aligner.align(ref, var)[0] # use first alignment only
    sub_list = [] # substitution (start, stop, ref, var)
    ins_list = [] # insertion
    del_list = [] # deletion
    pos = 1 # position on ref, 1-based
    #for r, v in zip(align.seqA, align.seqB):
    for r, v in zip(*alignment):
        if r == v:
            pos += 1
            continue
        elif r == '-': # insertion
            if ins_list == [] or pos > ins_list[-1][1]:
                ins_list.append([pos-1, pos, r, v, 'ins'])
            else: # continuous insertion
                #ins_list[-1][1] = pos
                ins_list[-1][2] += r
                ins_list[-1][3] += v
        elif v == '-': # deletion
            if del_list == [] or pos-1 > del_list[-1][1]:
                del_list.append([pos, pos, r, v, 'del'])
            else:
                del_list[-1][1] = pos
                del_list[-1][2] += r
                del_list[-1][3] += v
            pos += 1
        else: # substitution
            if sub_list == [] or pos-1 > sub_list[-1][1]:
                sub_list.append([pos, pos, r, v, 'sub'])
            else:
                sub_list[-1][1] = pos
                sub_list[-1][2] += r
                sub_list[-1][3] += v
            pos += 1
            
    mut_list = sorted(sub_list + ins_list + del_list, key=lambda x: x[0])
    n_var = len(mut_list)
    
    if n_var == 0:
        # no variants called
        mut_df = {
            'chr': [n_chr], 
            'start': [chr_pos], 
            'stop': [chr_pos+len(ref)-1], 
            'ref': [None], 
            'var': [None], 
            'type': ['wt'], 
            'HGVS_nomenclature': [None], 
            'var_seq': [var]
        }
        mut_df = pd.DataFrame(mut_df)
        return mut_df
    
    mut_list = list(zip(*mut_list))
    mut_df = {
        'chr': [n_chr] * n_var, 
        'start': mut_list[0], 
        'stop': mut_list[1], 
        'ref': mut_list[2], 
        'var': mut_list[3], 
        'type': mut_list[4]
    }
    mut_df = pd.DataFrame(mut_df)
    mut_df['start'] += chr_pos - 1
    mut_df['stop'] += chr_pos - 1
    
    # HGVS Nomenclature (https://varnomen.hgvs.org/)
    all_HGVS = []
    for _, row in mut_df.iterrows():
        start = row['start']
        stop = row['stop']
        r = row['ref']
        v = row['var']
        mut_type = row['type']
        
        if mut_type == 'sub':
            if len(v) == 1:
                # single nucleotide substitution
                HGVSG = '%s:g.%d%s>%s' % (n_chr, start, r, v)
            elif len(v) > 1:
                # multiple nucleotides substitution
                if r == revcomp(v):
                    # inversion
                    HGVSG = '%s:g.%d_%dinv' % (n_chr, start, stop)
                else:
                    # deletion-insertion
                    HGVSG = '%s:g.%d_%ddelins%s' % (n_chr, start, stop, v)
        elif mut_type == 'ins':
            start_0based = start - 1
            seq_ahead = ref[start_0based-len(v)+1 : start_0based+1]
            stop_0based = stop - 1
            seq_behind = ref[stop_0based : stop_0based+len(v)]
            if seq_ahead == v:
                # duplication (alignment condition 1)
                if len(v) == 1:
                    HGVSG = '%s:g.%ddup' % (n_chr, start)
                else:
                    HGVSG = '%s:g.%d_%ddup' % (n_chr, start-len(v)+1, start)
            elif seq_behind == v:
                # duplication (alignment condition 2)
                if len(v) == 1:
                    HGVSG = '%s:g.%ddup' % (n_chr, stop)
                else:
                    HGVSG = '%s:g.%d_%ddup' % (n_chr, stop, stop+len(v)-1)
            elif len(v) > 1 and seq_ahead == revcomp(v):
                # inverted duplication
                HGVSG = '%s:g.%d_%dins%d_%dinv' % (n_chr, start, stop, start-len(v)+1, start)
                #HGVSG = '%s:g.%d_%dins%d_%dinv' % (n_chr, stop+len(v)-1, stop+len(v), stop, stop+len(v)-1)
            else:
                # insertion
                HGVSG = '%s:g.%d_%dins%s' % (n_chr, start, stop, v)
        elif mut_type == 'del':
            HGVSG = '%s:g.%d_%ddel' % (n_chr, start, stop)
            
        all_HGVS.append(HGVSG)
    mut_df['HGVS_nomenclature'] = all_HGVS
    mut_df['var_seq'] = var
    return mut_df


def adp_trim_legacy(read_pair: tuple):
    """
    For ERBB2 179-plex and 226-plex panel (no rP spacer)
    Read structure: 
        on-target or short non-specific or dimer: 
            r1: spacer + UMI + amplicon + adapter
            r2: amplicon + UMI + spacer + adapter
        long non-specific:
            r1: spacer + UMI + amplicon (part)
            r2: spacer + amplicon (part)
    """
    r1, r2 = read_pair
    r1 = r1.decode('ascii')
    r2 = r2.decode('ascii')
    r1 = r1.split('\n')
    r2 = r2.split('\n')
    
    # line1 (record name)
    l1 = r1[0]
    name = l1[:l1.find(' ')] # content after space won't be kept in sam file

    #---------- TRIM CODE ----------#
    # line2 (sequence)
    l1 = r1[1]
    l2 = r2[1]
    
    is_error = 0 # 1: r1 & r2 merged
    is_matched = 1 # 1: r1 & r2 sequence exactly the same
    is_long = 1 # 1: longer than threshold after merging
    
    # merge r1 and r2
    # 1. find the head of r1 in r2
    pos_head = l2.find(revcomp(l1[SPACER_LEN:SPACER_SEARCH_LEN]))
    if pos_head != -1:
        end = SEARCH_LEN + pos_head + SPACER_LEN # end of amplicon in r1
        # NOTE: end may exceed max string index
        seq = l1[SPACER_LEN:end] # UMI + amplicon in r1
        length = end - SPACER_LEN
        exceed_len = length - len(seq)
        if exceed_len > 0:
            # make up with r2
            seq += revcomp(l2[:exceed_len])
        
        if seq != revcomp(l2[:end-SEARCH_LEN]):
            is_matched = 0
        
        name += '_UMI_%s' % (seq[:UMI_LEN])
        length -= UMI_LEN
        seq = seq[UMI_LEN:]
    else:
        # 2. find the tail of r1 in r2
        pos_tail = l2.find(revcomp(l1[-SEARCH_LEN:]))
        if pos_tail != -1: 
            seq = l1[UMI_SPACER_LEN:] + revcomp(l2[:pos_tail])
            length = len(seq)
            name += '_UMI_%s' % (l1[SPACER_LEN:UMI_SPACER_LEN])
        else:
            # cannot merge r1 & r2
            is_error = 1
    
    # length threshold
    if (not is_error) and (length < MIN_AMP_LEN_LEGACY):
        is_long = 0
    #---------- TRIM CODE ----------#
    
    # line 3 & 4
    if is_error:
        return None
    else:
        # line3
        line3 = r1[2]
        
        # line4 (quality)
        l1 = r1[3]
        l2 = r2[3]
        if pos_head != -1:
            qua = l1[UMI_SPACER_LEN:end]
            if exceed_len > 0:
                # make up with r2
                qua += l2[:exceed_len][::-1]
        elif pos_tail != -1:
            qua = l1[UMI_SPACER_LEN:] + l2[:pos_tail][::-1]
        else:
            pass
        
        outstr = '@%s\n%s\n%s\n%s\n' % (name, seq, line3, qua)
        return (
            outstr.encode('ascii'), 
            length, 
            is_long, 
            is_matched
        )


def adp_trim_short(read_pair: tuple):
    """
    Read structure: 
        on-target or short non-specific or dimer: 
            r1: spacer + UMI + amplicon + spacer + adapter
            r2: spacer + amplicon + UMI + spacer + adapter
        long non-specific:
            r1: spacer + UMI + amplicon (part)
            r2: spacer + amplicon (part)
    """
    r1, r2 = read_pair
    r1 = r1.decode('ascii')
    r2 = r2.decode('ascii')
    r1 = r1.split('\n')
    r2 = r2.split('\n')
    
    # line1 (record name)
    l1 = r1[0]
    name = l1[:l1.find(' ')] # content after space won't be kept in sam file

    #---------- TRIM CODE ----------#
    # line2 (sequence)
    l1 = r1[1]
    l2 = r2[1]
    
    is_error = 0 # 1: r1 & r2 merged
    is_matched = 1 # 1: r1 & r2 sequence exactly the same
    is_long = 1 # 1: longer than threshold after merging
    
    # merge r1 and r2
    # 1. find the head of r1 in r2
    pos_head = l2.find(revcomp(l1[SPACER_LEN:SPACER_SEARCH_LEN]))
    if pos_head != -1:
        end = SEARCH_LEN + pos_head # end of amplicon in r1
        # NOTE: end may exceed max string index
        seq = l1[SPACER_LEN:end] # UMI + amplicon in r1
        length = end - SPACER_LEN
        exceed_len = length - len(seq)
        if exceed_len > 0:
            # make up with r2
            seq += revcomp(l2[SPACER_LEN:SPACER_LEN+exceed_len])
        
        if seq != revcomp(l2[SPACER_LEN:end]):
            is_matched = 0
        
        name += '_UMI_%s' % (seq[:UMI_LEN])
        length -= UMI_LEN
        seq = seq[UMI_LEN:]
    else:
        # 2. find the tail of r1 in r2
        pos_tail = l2.find(revcomp(l1[-SEARCH_LEN:]))
        if pos_tail != -1: 
            seq = l1[UMI_SPACER_LEN:] + revcomp(l2[SPACER_LEN:pos_tail])
            length = len(seq)
            name += '_UMI_%s' % (l1[SPACER_LEN:UMI_SPACER_LEN])
        else:
            # cannot merge r1 & r2
            is_error = 1
    
    # length threshold
    if (not is_error) and (length < MIN_AMP_LEN):
        is_long = 0
    #---------- TRIM CODE ----------#
    
    # line 3 & 4
    if is_error:
        return None
    else:
        # line3
        line3 = r1[2]
        
        # line4 (quality)
        l1 = r1[3]
        l2 = r2[3]
        if pos_head != -1:
            qua = l1[UMI_SPACER_LEN:end]
            if exceed_len > 0:
                # make up with r2
                qua += l2[SPACER_LEN:SPACER_LEN+exceed_len][::-1]
        elif pos_tail != -1:
            qua = l1[UMI_SPACER_LEN:] + l2[SPACER_LEN:pos_tail][::-1]
        else:
            pass
        
        outstr = '@%s\n%s\n%s\n%s\n' % (name, seq, line3, qua)
        return (
            outstr.encode('ascii'), 
            length, 
            is_long, 
            is_matched
        )


def adp_trim_long(read_pair: tuple):
    """
    Read structure: 
        on-target or long non-specific:
            r1: spacer + UMI + amplicon (part)
            r2: spacer + amplicon (part)
        short non-specific or dimer: 
            r1: spacer + UMI + amplicon + spacer + adapter
            r2: spacer + amplicon + UMI + spacer + adapter
    """
    r1, r2 = read_pair
    r1 = r1.decode('ascii')
    r2 = r2.decode('ascii')
    r1 = r1.split('\n')
    r2 = r2.split('\n')
    
    # line1 (record name)
    l1 = r1[0]
    name = l1[:l1.find(' ')] # content after space won't be kept in sam file

    #---------- TRIM CODE ----------#
    # line2 (sequence)
    l1 = r1[1]
    l2 = r2[1]
    
    is_error = 0 # 1: r1 & r2 merged
    is_long = 1 # 1: longer than threshold after merging
    
    # merge r1 and r2
    # find the tail of r1 in r2
    pos_tail = l2.find(revcomp(l1[-SEARCH_LEN:]))
    if pos_tail != -1: 
        seq = l1[UMI_SPACER_LEN:] + revcomp(l2[SPACER_LEN:pos_tail])
        length = len(seq)
        name += '_UMI_%s' % (l1[SPACER_LEN:UMI_SPACER_LEN])
    
    else:
        # find the head of r1 in r2
        pos_head = l2.find(revcomp(l1[UMI_SPACER_LEN:UMI_SPACER_SEARCH_LEN]))
        if pos_head != -1:
            end = UMI_SEARCH_LEN + pos_head
            # NOTE: end may exceed max string index
            seq = l1[UMI_SPACER_LEN:end]
            length = end - UMI_SPACER_LEN
            exceed_len = length - len(seq)
            if exceed_len > 0:
                # make up with r2
                seq += revcomp(l2[SPACER_LEN:SPACER_LEN+exceed_len])
            name += '_UMI_%s' % (l1[SPACER_LEN:UMI_SPACER_LEN])
        
        else:
            is_error = 1

    # length threshold
    if (not is_error) and (length < MIN_AMP_LEN):
        is_long = 0
    #---------- TRIM CODE ----------#
    
    # line 3 & 4
    if is_error:
        return None
    else:
        # line3
        line3 = r1[2]
        
        # line4 (quality)
        l1 = r1[3]
        l2 = r2[3]

        if pos_tail != -1:
            qua = l1[UMI_SPACER_LEN:] + l2[SPACER_LEN:pos_tail][::-1]
        elif pos_head != -1:
            qua = l1[UMI_SPACER_LEN:end]
            if exceed_len > 0:
                # make up with r2
                qua += l2[SPACER_LEN:SPACER_LEN+exceed_len][::-1]
        else:
            pass
        
        outstr = '@%s\n%s\n%s\n%s\n' % (name, seq, line3, qua)
        return (
            outstr.encode('ascii'), 
            length, 
            is_long, 
            1
        )


def adp_trim_wrapper(r1_path: str, r2_path: str, out_path: str, dimer_path: str, sequencer: str, 
    panel_type: str, lib='lib', verbose=True, to_plot=True, plot_folder='', n_cpu=1, TRIM_BATCH_SIZE=2**30):
    """
    Wrapper for adp_trim_short and adp_trim_long
    :param r1_path: read 1 path, suffix should be .fastq.gz
    :param r2_path: read 2 path, suffix should be .fastq.gz
    :param out_path: output path, suffix is .gz.fastq
    :param dimer_path: output path for dimer reads, suffix is .gz.fastq
    :param sequencer: sequencing platform ('MiSeq' or 'else')
    :param panel_type: 'short' for short amplicon panel; 'long' for long amplicon panel
    :param lib: library name, for printing and plotting
    :param verbose: whether to print
    :param to_plot: whether to plot
    :param plot_folder: folder to save plots; default '' (save in current folder)
    :param n_cpu: number of CPUs for multiprocessing
    :param TRIM_BATCH_SIZE: number of bytes per batch; 2**31 takes ~10GB RAM (must be divisible by 4)
    :return: N, N_NoError, N_short (total reads, search succeeded reads, short reads)
    """
    if sequencer.lower() == 'miseq':
        read_delimiter = b'\n@M'
    elif sequencer.lower() == 'else':
        read_delimiter = b'\n@'
    else:
        raise ValueError('Invalid sequencer, must be "MiSeq" or "else".')

    if panel_type == 'short':
        trim_func = adp_trim_short
    elif panel_type == 'long':
        trim_func = adp_trim_long
    elif panel_type == 'legacy':
        trim_func = adp_trim_legacy
    else:
        raise ValueError('Invalid panel type. ')
    
    amp_len = [] # amplicon length only if r1 and r2 are successfully merged
    N_long = 0 # number of long amplicons (longer than MIN_AMP_LEN or MIN_AMP_LEN_LEGACY)
    N_matched = 0
    N_unmerged = 0
    N = 0 # total reads
    temp_r1 = b''
    temp_r2 = b''
    # gzip open with encoding will take longer time
    with gzip.open(r1_path, 'rb') as in_r1, gzip.open(r2_path, 'rb') as in_r2, \
    gzip.open(out_path, 'wb', compresslevel=1) as out, gzip.open(dimer_path, 'wb', compresslevel=1) as out_dimer:
        # r1 and r2 should have exactly the same number of bytes
        with multiprocessing.Pool(processes=n_cpu) as adp_trimmer:
            # get rid of thr first '@'
            in_r1.read(1)
            in_r2.read(1)

            tik = time()
            # read first batch
            r1_list = in_r1.read(TRIM_BATCH_SIZE)
            r2_list = in_r2.read(TRIM_BATCH_SIZE)

            while r1_list and r2_list:
                if r1_list[-1] == ord('\n'):
                    # avoid stopping between '\n' and '@'
                    r1_list += in_r1.read(1)
                    r2_list += in_r2.read(1)
                
                r1_list = r1_list.split(read_delimiter)
                r1_list[0] = temp_r1 + r1_list[0]
                temp_r1 = r1_list.pop(-1)
                
                r2_list = r2_list.split(read_delimiter)
                r2_list[0] = temp_r2 + r2_list[0]
                temp_r2 = r2_list.pop(-1)
                
                curr_reads_num = len(r1_list)
                if curr_reads_num != len(r2_list):
                    raise ValueError('read 1 and read 2 do not match!')
                N += curr_reads_num

                # trim current batch
                print('trimming current batch...')
                trimmed_reads = adp_trimmer.map(trim_func, zip(r1_list, r2_list))
                r1_list.clear()
                r2_list.clear()

                # go through trimmed reads
                # adding to '' or b'' could be very slow 
                long_lines = []
                dimer_lines = []
                for r in trimmed_reads:
                    if not r:
                        N_unmerged += 1
                        continue
                    amp_len.append(r[1]) # save amplicon length
                    N_long += r[2]
                    N_matched += r[3]
                    if r[2]:
                        # long amplicon
                        long_lines.append(r[0])
                    else:
                        # short amplicon
                        dimer_lines.append(r[0])
                trimmed_reads.clear()

                # dump to file
                print('dumping to file...')
                out.writelines(long_lines)
                out_dimer.writelines(dimer_lines)
                long_lines.clear()
                dimer_lines.clear()

                print('processed reads: %d, time passed: %.2fs' % (N, time()-tik))
                
                # read next batch
                r1_list = in_r1.read(TRIM_BATCH_SIZE)
                r2_list = in_r2.read(TRIM_BATCH_SIZE)

        # deal with the last read pair
        N += 1
        r = trim_func((temp_r1, temp_r2))
        if r:
            amp_len.append(r[1]) # save amplicon length
            N_long += r[2]
            N_matched += r[3]
            if r[2]:
                # long amplicon
                out.write(r[0])
            else:
                # short amplicon
                out_dimer.write(r[0])
        else:
            N_unmerged += 1
    
    N_NoError = len(amp_len)
    N_short = N_NoError - N_long
    N_unmatched = N_NoError - N_matched

    if verbose:
        print('total reads: %d' % N)
        print('unmerged: %d' % N_unmerged)
        print('unmatched: %d' % N_unmatched)
        print('search succeeded: %d, %.3f%%' % (N_NoError, 100*N_NoError/N))
        print('passed length threshold: %d, %.3f%%' % (N_NoError-N_short, 100*(N_NoError-N_short)/N_NoError))
        print('after trim rate: %.3f%%' % (100*(N_NoError-N_short)/N))
        
    plt.hist(amp_len, max(amp_len)-min(amp_len))
    if panel_type == 'legacy':
        plt.axvline(MIN_AMP_LEN_LEGACY, color='r', linewidth=1, linestyle='--')
    else:
        plt.axvline(MIN_AMP_LEN, color='r', linewidth=1, linestyle='--')
    plt.title(lib, size=14)
    plt.xlabel('amplicon length (bp)', size=14)
    plt.ylabel('# reads', size=14)
    savefig_path = os.path.join(plot_folder, '%s_length.png' % lib)
    plt.savefig(savefig_path, bbox_inches='tight')
    if to_plot:
        plt.show()
    else:
        plt.close()

    return N, N_NoError, N_short, N_unmerged, N_unmatched


def build_primer_hash(all_fp: list, all_rp: list):
    """
    Build up hash table for accelerated pairwise analysis
    :param all_fp: list of all forward primers
    :param all_rp: list of all reverse primers
    :return: hash table in dictionary format
    """
    HEAD_LEN = 5
    
    fp_dict = {}
    rp_dict = {}
    for i, (fp, rp) in enumerate(zip(all_fp, all_rp)):
        fp = fp.upper()
        rp = rp.upper()
        fp_key = fp[:HEAD_LEN]
        rp_key = rp[:HEAD_LEN]
        if not fp_key in fp_dict:
            fp_dict[fp_key] = [(i, fp, revcomp(fp))]
        else:
            fp_dict[fp_key].append((i, fp, revcomp(fp)))
        if not rp_key in rp_dict:
            rp_dict[rp_key] = [(i, rp, revcomp(rp))]
        else:
            rp_dict[rp_key].append((i, rp, revcomp(rp)))
        
    return {'HEAD_LEN': HEAD_LEN, 'tot_pair': i+1, 'fP': fp_dict, 'rP': rp_dict}


def primer_pairwise(primerhash: dict, amp_counter: dict, lib='lib'):
    """
    Analyze which pair of primer produced each amplicon
    :param primerhash: output of function build_primer_hash()
    :param amp_counter: a dictionary with amplicon sequence as key and its counts as value (e.g. Counter())
    :return: dataframe of each amplicon, matrix of count, dataframe of each primer pair
    """
    HEAD_LEN = primerhash['HEAD_LEN']
    tot_pair = primerhash['tot_pair']
    fp_dict = primerhash['fP']
    rp_dict = primerhash['rP']

    pair_cnt = Counter() # count of each primer pair
    # primer alignment (use default scores)
    aligner = PairwiseAligner(mode='local', match_score=1, mismatch_score=0, open_gap_score=0, extend_gap_score=0)
    all_pair = [] # list of primer index pair for each amplicon
    is_dimer = [] # whether this amplicon is dimer
    for seq, count in amp_counter.items():        
        # search the fP and rP that produced this amplicon
        head = seq[:HEAD_LEN]
        tail = revcomp(seq[-HEAD_LEN:])
        if (head in fp_dict) and (tail in rp_dict):
            fp_cand_list = fp_dict[head] # list of fP candidates
            rp_cand_list = rp_dict[tail] # list of rP candidates
            
            # find the best fP candidate
            if len(fp_cand_list) == 1:
                # only one primer with this heading, get its index
                i_fp, fp, _ = fp_cand_list[0]
            else:
                # more than one primers with this heading
                # get the local alignment score between each candidate primer and the amplicon
                # aln_score = [(pairwise2.align.localxx(seq, fp_cand, score_only=True), i_fp_cand, fp_cand) \
                #              for i_fp_cand, fp_cand, _ in fp_cand_list]
                aln_score = [(aligner.score(seq, fp_cand), i_fp_cand, fp_cand) \
                             for i_fp_cand, fp_cand, _ in fp_cand_list]
                # get the index of the primer with the best alignment score
                _, i_fp, fp = max(aln_score, key=lambda x: x[0])
                
            # find the best rP candidate
            if len(rp_cand_list) == 1:
                # only one primer with this heading, get its index
                i_rp, _, rp_rc = rp_cand_list[0]
            else:
                # more than one primers with this heading
                # get the local alignment score between each candidate primer and the amplicon
                # aln_score = [(pairwise2.align.localxx(seq, rp_cand_rc, score_only=True), i_rp_cand, rp_cand_rc) \
                #              for i_rp_cand, _, rp_cand_rc in rp_cand_list]
                aln_score = [(aligner.score(seq, rp_cand_rc), i_rp_cand, rp_cand_rc) \
                             for i_rp_cand, _, rp_cand_rc in rp_cand_list]
                # get the index of the primer with the best alignment score
                _, i_rp, rp_rc = max(aln_score, key=lambda x: x[0])
                
            pair_cnt[(i_fp, i_rp)] += count
            all_pair.append((i_fp, i_rp))

            # identify dimers based on length
            if len(fp) + len(rp_rc) >= len(seq) + 3:
                is_dimer.append(True)
            else:
                is_dimer.append(False)
        
        else:
            all_pair.append((None, None))
            is_dimer.append(None)

    # matrix of counts
    mtx_count = np.zeros([tot_pair, tot_pair])
    for pair, count in pair_cnt.items():
        mtx_count[pair] = count

    # dataframe of each amplicon (sequence, count, rario, fp index, rp index)
    N = sum(amp_counter.values())
    fp_list, rp_list = list(zip(*all_pair))
    df_amp = pd.DataFrame({
            'seq': amp_counter.keys(), 
            'count': amp_counter.values(), 
            'ratio': np.array(list(amp_counter.values()))/N, 
            'fP': fp_list, 
            'rP': rp_list, 
            'is_dimer': is_dimer
        })
    df_amp.sort_values('count', axis='index', ascending=False, inplace=True)
    df_amp.reset_index(drop=True, inplace=True)

    # dataframe of each primer pair (fp index, rp index, count, rario)
    matched = sum(pair_cnt.values())
    fp_list, rp_list = list(zip(*pair_cnt.keys()))
    df_pair = pd.DataFrame({
            'fP': fp_list, 
            'rP': rp_list, 
            'count': pair_cnt.values(), 
            'ratio': np.array(list(pair_cnt.values()))/matched
        })
    df_pair.sort_values('count', axis='index', ascending=False, inplace=True)
    df_pair.reset_index(drop=True, inplace=True)

    return df_amp, mtx_count, df_pair


def UMI_analysis(args):
    """
    UMI analysis for single amplicon
    :param sorted_bam_path: path to a indexed, sorted bam file
    :param ref_row: row of the design sheet, corresponding to the amplicon
    :param raw_path: path to save raw UMI grouping
    :param plot_path: path to save family size distribution plots
    :param voted_path: path to save the voted sequence for each UMI family
    :param vaf_path: path to save the VAF of each unique voted sequence
    :param save_raw_UMI: True to save raw UMI grouping file
    :param DS_rate: downsampling rate
    :return: tot_umi (total UMI count), ref_mut_df (all variants called), N_GroupingLoss, N_VotingLoss
    """
    _UMI_FILTER = 'dynamic' # 'dynamic' or 'gaussian', for internal use only

    BIN_WIDTH = 0.35
    GAUSSIAN_CUTOFF = 2
    GAUSS_SIGMA = 1
    #p0 = [300., 5., 1.]
    N_SIGMA = 3

    sorted_bam_path, ref_row, raw_path, plot_path, voted_path, vaf_path, save_raw_UMI, DS_rate = args

    samfile = pysam.AlignmentFile(sorted_bam_path, 'rb')

    ref_name = ref_row['plex_id']
    ref_chr = ref_row['chr']
    ref_pos = ref_row['insert_start']
    ref_seq = ref_row['insert_in'].upper()
    fp_val_seq = ref_row['fP'][-PRIMER_VAL_LEN:].upper()
    rp_val_seq = revcomp(ref_row['rPin'][-PRIMER_VAL_LEN:]).upper()
    fp_len = len(ref_row['fP'])
    rp_len = len(ref_row['rPin'])
    
    #---------- Group Reads by UMI ----------#
    # group sequences (reads) with the same UMI
    N_GroupingLoss = 0 # reads loss during grouping
    umi_dict = dict()
    for record in samfile.fetch(contig=ref_name):
        if np.random.rand() > DS_rate:
            continue

        # primer validation
        if record.query_sequence[fp_len-PRIMER_VAL_LEN:fp_len] != fp_val_seq or record.query_sequence[-rp_len:-rp_len+PRIMER_VAL_LEN] != rp_val_seq:
            # cannot locate primer
            N_GroupingLoss += 1
            continue

        # save insert sequence and quality
        umi = record.query_name[-UMI_LEN:]
        seq = record.query_sequence[fp_len:-rp_len]
        if not seq:
            # in case sequence is empty
            N_GroupingLoss += 1
            continue
        qual = tuple(record.query_qualities)[fp_len:-rp_len]
        if umi in umi_dict:
            umi_dict[umi].append((seq, qual))
        else:
            umi_dict[umi] = [(seq, qual)]

    # save umi_dict to json
    if DS_rate >= 1:
        json_path = os.path.join(raw_path, '%s.json' % ref_name)
        if save_raw_UMI:
            with open(json_path, 'w') as f:
                json.dump(umi_dict, f, indent=4)
    #---------- Group Reads by UMI ----------#


    #---------- UMI Filtering ----------#
    # number of raw UMI families
    N_UMI_raw = len(umi_dict)

    # first filtering: remove UMIs containing 'G'
    umi_to_del = [umi for umi in umi_dict if 'G' in umi]
    for umi in umi_to_del:
        del umi_dict[umi]
    N_UMI_NoG = len(umi_dict)
    
    # second filtering: dynamic cutoff
    family_size = [len(family) for family in umi_dict.values()] # list of family sizes
    # determine family size threshold
    if (_UMI_FILTER == 'dynamic') and (len(umi_dict) >= TOP_N):
        topn_count = sorted(family_size, reverse=True)[:TOP_N]
        dynamic_th = math.floor(DYNAMIC_RATE * sum(topn_count)/TOP_N)
        final_th = max(STATIC_TH, dynamic_th)
    elif (_UMI_FILTER == 'gaussian') and family_size:
        family_size_log = np.log2(family_size)
        # number of bins must be greater than 1
        hist_height, bin_edge = np.histogram(family_size_log, bins=int(max(family_size_log)/BIN_WIDTH)+1)
        bin_center = (bin_edge[:-1] + bin_edge[1:])/2
        hist_filter = bin_center > GAUSSIAN_CUTOFF # only fit part of data (do not include UMI families with small sizes)
        bin_center_filtered = bin_center[hist_filter]
        hist_height_filtered = hist_height[hist_filter]
        # determine family size threshold by fitting to gaussian distribution
        if len(bin_center_filtered) >= 3: # cannot converge if only one bin
            try:
                gauss_idx = np.argmax(hist_height_filtered)
                p0 = [hist_height_filtered[gauss_idx], bin_center_filtered[gauss_idx], GAUSS_SIGMA]
                coeff, _ = curve_fit(gauss, bin_center_filtered, hist_height_filtered, p0=p0)
                dynamic_th = 2**(coeff[1] - N_SIGMA*abs(coeff[2]))
                final_th = max(STATIC_TH, dynamic_th)
                hist_height_fit = gauss(bin_center_filtered, *coeff)
            except RuntimeError:
                # there still could be cases where the fitting cannot converge
                final_th = STATIC_TH
        else:
            final_th = STATIC_TH
    else:
        final_th = STATIC_TH
    # remove UMIs with family size under threshold
    umi_to_del = [umi for umi, family in umi_dict.items() if len(family) <= final_th]
    for umi in umi_to_del:
        del umi_dict[umi]
    
    # total number of UMIs after filtering
    tot_umi = len(umi_dict)

    if DS_rate < 1:
        umi_dict.clear()
        return tot_umi

    # plot family size distribution
    savefig_path = os.path.join(plot_path, '%s.svg' % ref_name)
    fig = plt.figure(figsize=(6, 5))
    if (_UMI_FILTER == 'dynamic') and family_size:
        family_size_log = np.log2(family_size)
        if max(family_size) != min(family_size):
            #plt.hist(family_size_log, max(family_size)-min(family_size), log=False)
            hist_arr = plt.hist(family_size_log, int(max(family_size_log)/0.35), log=False)
            plt.ylim(top=1.1*max(hist_arr[0][1:])) # ignore UMIs with family size = 1
        else:
            plt.hist(family_size_log, log=False)
    elif (_UMI_FILTER == 'gaussian') and family_size:
        plt.hist(family_size_log, bins=bin_edge, log=False)
        if len(hist_height) > 1:
            plt.ylim(top=1.1*max(hist_height[1:])) # ignore the scale of UMIs with family size = 1
        try:
            plt.plot(bin_center_filtered, hist_height_fit, color='orange', linewidth=1)
            plt.axvline(coeff[1], color='k', linewidth=1, linestyle='--')
            plt.axvline(coeff[1] - abs(coeff[2]), color='k', linewidth=1, linestyle='--')
            plt.axvline(coeff[1] - 2*abs(coeff[2]), color='k', linewidth=1, linestyle='--')
            plt.axvline(coeff[1] - 3*abs(coeff[2]), color='k', linewidth=1, linestyle='--')
        except NameError:
            pass
    
    if family_size:
        plt.axvline(np.log2(final_th), color='r', linewidth=1, linestyle='--')
        plt.title(ref_name, fontsize=14)
        text_str = 'med=%d\nmax=%d\nraw:%d\nno G:%d\nfinal:%d' % \
        (np.median(family_size), max(family_size), N_UMI_raw, N_UMI_NoG, tot_umi)
        plt.annotate(text_str, (0.65, 0.7), xycoords='figure fraction', fontsize=14)
        plt.xlabel('log2(UMI family size)', fontsize=14)
        plt.ylabel('# UMI family', fontsize=14)
    else: # no UMI before filtering
        plt.annotate('No UMI found', (0.35, 0.5), xycoords='figure fraction', fontsize=14)
        plt.xticks([])
        plt.yticks([])
        plt.title(ref_name, fontsize=14)
    plt.savefig(savefig_path, bbox_inches='tight')
    plt.close()
    #---------- UMI Filtering ----------#


    #---------- UMI Voting & VAF ----------#
    # single base voting, weighted by sequencing quality
    # collapse sequences in the same UMI family into one voted sequence
    N_VotingLoss = 0 # reads loss during voting
    voted_dict = dict()
    for umi, family in umi_dict.items():
        # most common sequence length in this UMI family
        most_common_len = Counter([len(s) for s, _ in family]).most_common(1)[0][0]
        # select family member(s) with the same length; convert to complex format
        family_same_len = [(seq2complex(s), q) for s, q in family if len(s) == most_common_len]
        N_VotingLoss += len(family) - len(family_same_len)
        # sequence matrix & quality matrix
        seq_mtx, qual_mtx = list(zip(*family_same_len))
        seq_mtx = np.array(seq_mtx)
        qual_mtx = np.array(qual_mtx)
        # weighted by quality
        seq_voted_complex = np.sum(np.multiply(seq_mtx, qual_mtx), axis=0)
        # convert back to string format
        seq_voted = complex2seq(seq_voted_complex)
        voted_dict[umi] = seq_voted

    # save voted_dict to json
    json_path = os.path.join(voted_path, '%s.json' % ref_name)
    with open(json_path, 'w') as f:
        json.dump(voted_dict, f, indent=4)
    
    # calculate VAF
    # collapse UMIs with the same voted sequence, to calculate VAF
    vaf_dict = Counter(voted_dict.values()).most_common()
    vaf_dict = {seq: [count, count/tot_umi] for seq, count in vaf_dict}

    # save vaf_dict to json
    json_path = os.path.join(vaf_path, '%s.json' % ref_name)
    with open(json_path, 'w') as f:
        json.dump(vaf_dict, f, indent=4)
    #---------- UMI Voting & VAF ----------#


    #---------- Variant Call ----------#
    # some amplicons are designed on the minus strand for 179-plex and 226-plex panel
    try:
        strand = ref_row['strand']
        if strand == 2:
            ref_pos = ref_row['insert_stop']
            ref_seq = revcomp(ref_seq)
    except KeyError:
        strand = 1

    ref_mut_df = []
    for i, (var_seq, (mol_count, vaf)) in enumerate(vaf_dict.items()):
        if strand == 2:
            var_seq = revcomp(var_seq)
        mut_df = variant_call(ref_seq, var_seq, ref_chr, ref_pos)
        mut_df['mol_group'] = i+1
        mut_df['mol_count'] = mol_count
        mut_df['VAF'] = vaf
        ref_mut_df.append(mut_df)
    if ref_mut_df != []: # must use is when comparing against []
        ref_mut_df = pd.concat(ref_mut_df, axis='index', ignore_index=True)
        ref_mut_df['plex_id'] = ref_name
    else:
        # current amplicon (ref_name) has no UMI
        ref_mut_df = None # must set to None instead of leaving it as []
        # since pd.DataFrame cannot compare with [] (Unable to coerce to Series)
    #---------- Variant Call ----------#

    # flush memory (use clear() instead of del)
    umi_dict.clear()
    voted_dict.clear()
    vaf_dict.clear()
    
    return tot_umi, ref_mut_df, N_GroupingLoss, N_VotingLoss


def UMI_analysis_wrapper(sorted_bam_path, df_design, hp_dict, lib, umi_folder, DS_folder, mutation_folder, n_cpu=1, save_raw_UMI=True):
    """
    Wrapper function for UMI_analysis
    :param sorted_bam_path: path to a indexed, sorted bam file
    :param df_design: pd.DataFrame of the design sheet
    :param hp_dict: dict of homopolymer coordinates
    :param lib: name of the library
    :param umi_folder: path to save intermediate UMI analysis files
    :param mutation_folder: path to save variants called
    :param n_cpu: number of core to do multiprocessing
    :param save_raw_UMI: True to save raw UMI grouping file
    :return: umi_num (dict of UMI count of each amplicon), N_GroupingLoss, N_VotingLoss
    """
    df_design = list(df_design.iterrows())
    df_design = list(zip(*df_design))[1] # list of rows
    N_ref = len(df_design)

    # create output folders
    raw_path = os.path.join(umi_folder, '%s_UMI_1-raw' % lib)
    plot_path = os.path.join(umi_folder, '%s_FamilySize' % lib)
    voted_path = os.path.join(umi_folder, '%s_UMI_2-voted' % lib)
    vaf_path = os.path.join(umi_folder, '%s_UMI_3-VAF' % lib)
    os.system('mkdir %s' % raw_path)
    os.system('mkdir %s' % plot_path)
    os.system('mkdir %s' % voted_path)
    os.system('mkdir %s' % vaf_path)

    # analysis UMI across multiple cores
    with multiprocessing.Pool(processes=n_cpu) as analyzer:
        func_input = zip(
            [sorted_bam_path]*N_ref, 
            df_design, 
            [raw_path]*N_ref, 
            [plot_path]*N_ref, 
            [voted_path]*N_ref, 
            [vaf_path]*N_ref, 
            [save_raw_UMI]*N_ref, 
            [1]*N_ref # downsampling rate = 1 (no dowmsampling)
        )
        umi_results = analyzer.map(UMI_analysis, func_input)
        umi_num, df_AllVar, N_GroupingLoss, N_VotingLoss = list(zip(*umi_results))

        # calculate CV of downsampling
        print('downsampling...')
        DS_mtx = [umi_num]
        for DS_rate in DS_RATE_LS:
            func_input = zip(
                [sorted_bam_path]*N_ref, 
                df_design, 
                [raw_path]*N_ref, 
                [plot_path]*N_ref, 
                [voted_path]*N_ref, 
                [vaf_path]*N_ref, 
                [save_raw_UMI]*N_ref, 
                [DS_rate]*N_ref
            )
            umi_num_DS = analyzer.map(UMI_analysis, func_input)
            DS_mtx.append(umi_num_DS)
    DS_mtx = np.array(DS_mtx)
    DS_mean = np.mean(DS_mtx, axis=0)
    DS_mean += np.array([int(n == 0) for n in DS_mean])
    DS_mask = 2 * (np.array([int(n > CV_MIN_UMI_TH) for n in DS_mean]) - 0.5) # amplicons with too few UMIs will have negative CV
    DS_std = np.std(DS_mtx, axis=0, ddof=1)
    DS_cv = DS_mask * DS_std / DS_mean
    CV_th = np.std(DS_RATE_LS, ddof=1) / np.mean(DS_RATE_LS) # CV threshold
    DS_cv_filter = [cv for cv in DS_cv if cv >= 0] # ignore amplicons with low UMI count
    N_low_CV = len([cv for cv in DS_cv_filter if cv < CV_th])

    # save downsampling UMI count
    np.savetxt(os.path.join(DS_folder, '%s_DS.csv' % lib), DS_mtx, delimiter=',')

    # plot downsampling CV
    savefig_path = os.path.join(DS_folder, '%s_CV.svg' % lib)
    fig = plt.figure(figsize=(6, 5))
    plt.hist(DS_cv_filter)
    plt.axvline(CV_th, color='r', linewidth=1, linestyle='--')
    plt.title(lib, fontsize=14)
    text_str = 'total=%d\nshown=%d\nlowCV:%d' % (len(DS_cv), len(DS_cv_filter), N_low_CV)
    plt.annotate(text_str, (0.65, 0.7), xycoords='figure fraction', fontsize=14)
    plt.xlabel('CV %s' % DS_RATE_LS, fontsize=14)
    plt.ylabel('# amplicon', fontsize=14)
    plt.savefig(savefig_path, bbox_inches='tight')
    plt.close()

    umi_num = {row['plex_id']: n for row, n in zip(df_design, umi_num)}
    DS_cv = {row['plex_id']: n for row, n in zip(df_design, DS_cv)}
    N_GroupingLoss = sum(N_GroupingLoss)
    N_VotingLoss = sum(N_VotingLoss)

    # remove empty dataframe
    df_AllVar = [df for df in df_AllVar if df is not None]
    if df_AllVar == []:
        # this library has no UMI at all
        return umi_num, N_GroupingLoss, N_VotingLoss, DS_cv, DS_mtx
    df_AllVar = pd.concat(df_AllVar, axis='index', ignore_index=True)

    # filter variants
    df_filter = df_AllVar[df_AllVar['var'] != 'N'] # drop 'N'
    df_filter = df_filter[df_filter['type'] != 'wt'] # drop WT
    df_filter = df_filter[df_filter['mol_count'] >= UMI_COUNT_TH] # remove low UMI counts
    # remove homopolymer indel
    df_indel = df_filter[df_filter['type'] != 'sub']
    hp_idx = [] # index of rows to drop
    for _, row in df_indel.iterrows():
        hp_list = hp_dict[row['plex_id']]
        hp_flag = False
        for hp_coords in hp_list:
            if hp_coords[0] <= row['start'] <= hp_coords[1]:
                hp_flag = True
                hp_idx.append(row.name)
                break
    df_hp = df_filter.loc[hp_idx, ] # record homopolymer indel
    df_filter.drop(index=hp_idx, inplace=True)

    # wirte variants to excel
    with pd.ExcelWriter('%s/%s_varlist.xlsx' % (mutation_folder, lib)) as writer:
        df_filter.to_excel(writer, sheet_name='final', index=True)
        df_hp.to_excel(writer, sheet_name='homopolymer', index=True)
        df_AllVar.to_excel(writer, sheet_name='raw', index=True)

    print('compressing intermediate files...')
    if save_raw_UMI:
        make_tarfile(raw_path, compress_type='gz', compresslevel=1, keep_source=False)
    else:
        # delete the empty folder
        os.system('rm -r %s' % raw_path)
    make_tarfile(plot_path, compress_type='gz', compresslevel=1, keep_source=False)
    make_tarfile(voted_path, compress_type='gz', compresslevel=1, keep_source=False)
    make_tarfile(vaf_path, compress_type='gz', compresslevel=1, keep_source=False)

    return umi_num, N_GroupingLoss, N_VotingLoss, DS_cv, DS_mtx


def build_index(design_file_path, n_tube):
    """
    generate index files, including homopolymers, bowtie2, primer hash
    :param design_file_path: path to the design spreadsheet
    :param n_tube: number of tubes
    """
    HP_LEN_TH = 4 # homopolymer length threshold

    # create index folder
    index_folder = 'index'
    os.system('mkdir %s' % index_folder)

    for i_tube in range(1, n_tube+1):
        df = pd.read_excel(design_file_path, 'tube_%d' % i_tube, index_col=0, engine='openpyxl')
        if list(df.index) != list(range(len(df))):
            raise ValueError('Invalid index column in design sheet tube_%d, index must be 0 to n.' % i_tube)
        
        # find homopolymer
        hp_dict = {}
        for _, row in df.iterrows():
            ref_name = row['plex_id']
            insert_seq = row['insert_in'].upper()
            insert_pos = row['insert_start']
            
            hp_list = [] # homopolymer coordinates of current insert
            b_before = insert_seq[0]
            chain_start = insert_pos
            chain_stop = insert_pos
            for i, b in enumerate(insert_seq[1:]):
                curr_pos = i + insert_pos + 1
                if b == b_before:
                    chain_stop = curr_pos
                else:
                    if chain_stop - chain_start + 1 >= HP_LEN_TH:
                        # homopolymer longer than threshold
                        hp_list.append([chain_start, chain_stop])
                    chain_start = curr_pos
                b_before = b
            if chain_stop - chain_start + 1 >= HP_LEN_TH:
                # handle the last homopolymer
                hp_list.append([chain_start, chain_stop])
            hp_dict[ref_name] = hp_list
        with open('%s/homopolymer_tube%d.json' % (index_folder, i_tube), 'w') as f:
            json.dump(hp_dict, f, indent=4)
        
        # bowtie2 index
        with open('%s/tube_%d.fasta' % (index_folder, i_tube), 'w') as f:
            for _, row in df.iterrows():
                f.write('>%s\n' % row['plex_id'])
                f.write('%s\n' % row['amp_in'].upper())
        os.system('bowtie2-build %s/tube_%d.fasta %s/tube_%d' % (index_folder, i_tube, index_folder, i_tube))
        
        # primer hash
        all_fp = list(df['fP'])
        all_rp = list(df['rPin'])
        primerhash = build_primer_hash(all_fp, all_rp)
        with open('%s/PrimerHeadDict_tube%d.json' % (index_folder, i_tube), 'w') as f:
            json.dump(primerhash, f, indent=4)
        
        # save a copy of the design spread sheet
        df.to_csv('%s/design_tube%d.csv' % (index_folder, i_tube), index=True)

    print('Index built successfully. ')


def run_pipeline(lib_list: list, sequencer: str, r1_path_format: str, r2_path_format: str, panel_type='long', 
    total_tube=2, verbose=False, to_plot=False, keep_sam=False, save_raw_UMI=True, performace_level=1, 
    TRIM_ALIGN_SWITCH=True, UMI_SWITCH=True, UNIFORMITY_SWITCH=True, DIMER_SWITCH=True, NS_SWITCH=True):
    """
    Pipeline to link all analysis together
    :param lib_list: list of library names
    :param r1_path_format: path to raw read 1 (.fastq.gz), formated string (e.g., 'raw_reads/%s_R1.fastq.gz')
    :param r1_path_format: path to raw read 2 (.fastq.gz), formated string (e.g., 'raw_reads/%s_R2.fastq.gz')
    :param panel_type: 'long' or 'short'
    :param total_tube: int, total number of tubes
    :param verbose: True to print outputs
    :param to_plot: True to show plots
    :param keep_sam: True to keep sam files
    :param save_raw_UMI: True to save raw UMI grouping file
    :param performace_level: int, 1 to use fewer RAM and half of all logic cores, 2 to use more RAM and all cores
    :param TRIM_ALIGN_SWITCH: True to do adapter trim and alignment
    :param UMI_SWITCH: True to do UMI analysis
    :param UNIFORMITY_SWITCH: True to output reads uniformity
    :param DIMER_SWITCH: True to do primer dimer analysis
    :param NS_SWITCH: True to do non-specific amplification analysis
    :return: 
    """
    if performace_level == 1:
        TRIM_BATCH_SIZE = 2**30
        n_cpu = multiprocessing.cpu_count()//2
    elif performace_level == 2:
        TRIM_BATCH_SIZE = 2**31
        n_cpu = multiprocessing.cpu_count()
    else:
        raise ValueError('Invalid performace_level')
    print('Performance level %d. Running on %d CPU(s)...\n' % (performace_level, n_cpu))

    # create folders
    index_folder = 'index'
    trim_folder = 'adapter_trim'
    bowtie2_folder = 'bowtie2_out'
    dimer_folder = 'dimer'
    sam_folder = 'sam_files'
    bam_folder = 'bam_files'
    umi_folder = 'UMI'
    DS_folder = 'downsampling'
    mutation_folder = 'mutation'
    uniformity_folder = 'uniformity'
    NS_folder = 'non-specific'
    if TRIM_ALIGN_SWITCH:
        os.system('mkdir %s' % trim_folder)
        os.system('mkdir %s' % bowtie2_folder)
        os.system('mkdir %s' % sam_folder)
        os.system('mkdir %s' % bam_folder)
        os.system('mkdir %s' % dimer_folder)
    if UMI_SWITCH:
        os.system('mkdir %s' % umi_folder)
        os.system('mkdir %s' % DS_folder)
        os.system('mkdir %s' % mutation_folder)
    if UNIFORMITY_SWITCH:
        os.system('mkdir %s' % uniformity_folder)
    if DIMER_SWITCH:
        pass
    if NS_SWITCH:
        os.system('mkdir %s' % NS_folder)

    # metrics
    all_total = [] # total reads of all libraries
    all_noerror = [] # search succeeded reads of all libraries
    all_short = [] # short reads of all libraries
    all_unmerged = []
    all_unmatched = []
    all_aligned = [] # aligned reads of all libraries
    all_mad = [] # Median absolute deviation
    all_max_min = []
    all_95_5 = []
    all_n_dropout = []
    all_med_UMI = []

    umi_num_dict = [{} for _ in range(total_tube)] # UMI number of all tubes
    DS_cv_dict = [{} for _ in range(total_tube)] # downsampling CV of all tubes
    read_dict = [{} for _ in range(total_tube)] # reads per plex of all tubes
    # analyze each library
    for i_lib, lib in enumerate(lib_list):
        print('----------------------------------------------')
        print(lib)
        tik = time()

        try:
            n_tube = int(lib[0])
        except ValueError:
            n_tube = 1
        r1_path = r1_path_format % lib
        r2_path = r2_path_format % lib
        trim_path = '%s/%s_trimmed.fastq.gz' % (trim_folder, lib)
        bowtie2_path = '%s/%s_bowtie2.txt' % (bowtie2_folder, lib) # bowtie2 output
        dimer_path = '%s/%s_dimer.fastq.gz' % (dimer_folder, lib)
        sam_path = '%s/%s.sam' % (sam_folder, lib)
        bam_path = '%s/%s_sorted.bam' % (bam_folder, lib)
        
        #------------------- adapter trim and alignment -------------------#
        if TRIM_ALIGN_SWITCH:
            print('adapter trim...')
            tik1 = time()
            with futures.ProcessPoolExecutor(max_workers=1) as executor:
                # create child process to prevent building up memory
                args = (r1_path, r2_path, trim_path, dimer_path, sequencer, 
                    panel_type, lib, verbose, to_plot, trim_folder, n_cpu, TRIM_BATCH_SIZE)
                N, N_NoError, N_short, N_unmerged, N_unmatched = executor.submit(adp_trim_wrapper, *args).result()
            print('adapter trim finished in %.2fs\n' % (time()-tik1))

            # align with bowtie2
            tik1 = time()
            print('aligning to tube %d reference...' % n_tube)
            cmd = '(bowtie2 -p %d -x %s/tube_%d -U %s -S %s) 2>%s' % (n_cpu, index_folder, n_tube, trim_path, sam_path, bowtie2_path)
            os.system(cmd)
            # parse bowtie2 results
            with open(bowtie2_path) as f:
                for line in f:
                    if 'aligned exactly 1 time' in line:
                        break
                #lines = f.readlines()
            aligned = int(line.strip().split(' ')[0])
            print('aligned: %d, %.3f%%' % (aligned, 100*aligned/(N_NoError-N_short)))
            print('on-target rate: %.3f%%' % (100*aligned/N))
            # sort and index sam file with pysam
            print('converting sam to bam...')
            pysam.sort('-@', str(n_cpu), '-o', bam_path, sam_path)
            pysam.index('-@', str(n_cpu), bam_path)
            print('alignment finished in %.2fs\n' % (time()-tik1))
            
            all_total.append(N)
            all_noerror.append(N_NoError)
            all_short.append(N_short)
            all_unmerged.append(N_unmerged)
            all_unmatched.append(N_unmatched)
            all_aligned.append(aligned)
        #------------------- adapter trim and alignment -------------------#
        
        
        #------------------- UMI and variant call -------------------#
        if UMI_SWITCH:
            # tik1 = time()
            # with futures.ProcessPoolExecutor(max_workers=1) as executor:
            #     # create child process to prevent building up memory
            #     args = (bam_path, lib, n_tube, index_folder, umi_folder, mutation_folder, save_raw_UMI)
            #     umi_num = executor.submit(UMI_analysis_wrapper, *args).result()
            # print('UMI analysis finished in %.2fs\n' % (time()-tik1))

            tik1 = time()
            # load panel sequences and coordinates
            df_design = pd.read_csv('%s/design_tube%d.csv' % (index_folder, n_tube), index_col=0)
            # load homopolymer coordinates
            with open('%s/homopolymer_tube%d.json' % (index_folder, n_tube)) as f:
                hp_dict = json.load(f)
            print('UMI analysis (UMI couting & variant call)...')
            umi_num, N_GroupingLoss, N_VotingLoss, DS_cv, DS_mtx = UMI_analysis_wrapper(bam_path, df_design, hp_dict, 
                lib, umi_folder, DS_folder, mutation_folder, n_cpu, save_raw_UMI)
            med_UMI_count = np.median(list(umi_num.values()))
            print('grouping loss: %d, voting loss: %d' % (N_GroupingLoss, N_VotingLoss))
            print('median UMI count: %d' % med_UMI_count)
            print('UMI analysis finished in %.2fs\n' % (time()-tik1))

            # save and plot UMI counts
            umi_num_dict[n_tube-1]['amplicon_name'] = list(umi_num.keys())
            umi_num_dict[n_tube-1][lib] = list(umi_num.values())
            all_med_UMI.append(med_UMI_count)
            umi_num_sorted = sorted(umi_num.values())
            plt.plot(umi_num_sorted)
            plt.axhline(6000, color='k', linewidth=1, linestyle='--')
            plt.axhline(np.median(umi_num_sorted), color='k', linewidth=1, linestyle='--')
            plt.xlabel('amplicons (sorted)', size=14)
            plt.ylabel('# UMIs', size=14)
            plt.title(lib, size=14)
            plt.savefig('%s/%s_UMI.png' % (umi_folder, lib), bbox_inches='tight')
            if to_plot:
                plt.show()
            else:
                plt.close()

            # save downsampling CV
            DS_cv_dict[n_tube-1]['amplicon_name'] = list(DS_cv.keys())
            DS_cv_dict[n_tube-1][lib] = list(DS_cv.values())

            #-------------------------------------#
            # added for 105-plex RNA panel
            # umi_num_sorted = umi_num.values()
            # plt.figure(figsize=[20, 5])
            # plt.bar(range(len(umi_num_sorted)), umi_num_sorted, 0.4)
            # plt.grid(axis='y')
            # for i in range(len(umi_num_sorted)//5 - 1):
            #     plt.axvline(5*(i+1)-0.5, color='k', linewidth=1, linestyle='--')
            # plt.xlabel('amplicons', size=14)
            # plt.ylabel('# UMIs', size=14)
            # plt.title(lib, size=14)
            # plt.savefig('%s/%s_UMI_unsorted.png' % (umi_folder, lib), bbox_inches='tight')
            # if to_plot:
            #     plt.show()
            # else:
            #     plt.close()

            # umi_num_sorted = np.log10(list(umi_num.values()))
            # plt.figure(figsize=[20, 5])
            # plt.bar(range(len(umi_num_sorted)), umi_num_sorted, 0.4)
            # plt.grid(axis='y')
            # for i in range(len(umi_num_sorted)//5 - 1):
            #     plt.axvline(5*(i+1)-0.5, color='k', linewidth=1, linestyle='--')
            # plt.xlabel('amplicons', size=14)
            # plt.ylabel('log10(# UMIs)', size=14)
            # plt.title(lib, size=14)
            # plt.savefig('%s/%s_UMI_unsorted_log10.png' % (umi_folder, lib), bbox_inches='tight')
            # if to_plot:
            #     plt.show()
            # else:
            #     plt.close()
            #-------------------------------------#
        else:
            all_med_UMI.append(-1)
            
        #------------------- UMI and variant call -------------------#
        
        
        #------------------- uniformity -------------------#
        if UNIFORMITY_SWITCH:
            samfile = pysam.AlignmentFile(bam_path, 'rb')
            stats = samfile.get_index_statistics()
            mapped = np.array([s.mapped for s in stats])
            amp_name = [s.contig for s in stats]
            
            read_dict[n_tube-1]['amplicon_name'] = amp_name
            read_dict[n_tube-1][lib] = mapped
            
            mapped_log10 = np.log10(mapped + 1) # aviod log10(0)
            mapped_log10_sorted = np.sort(mapped_log10)
            mapped_log10_dev = np.absolute(mapped_log10 - np.median(mapped_log10))
            all_mad.append(np.median(mapped_log10_dev))
            all_max_min.append((max(mapped_log10)-min(mapped_log10))/max(mapped_log10))
            i_95 = math.floor(0.95 * len(mapped_log10))
            i_5 = math.floor(0.05 * len(mapped_log10))
            #all_95_5.append((mapped_log10_sorted[i_95]-mapped_log10_sorted[i_5])/mapped_log10_sorted[i_95])
            all_95_5.append(mapped_log10_sorted[i_95]-mapped_log10_sorted[i_5])
            all_n_dropout.append(len([n for n in mapped_log10 if n == 0]))

            plt.plot(mapped_log10_sorted)
            plt.axhline(max(mapped_log10), color='k', linewidth=1, linestyle='--')
            plt.axhline(np.median(mapped_log10), color='k', linewidth=1, linestyle='--')
            plt.axhline(min(mapped_log10), color='k', linewidth=1, linestyle='--')
            #plt.ylim([0, max(mapped_log10)+1])
            plt.ylim([0, 6])
            plt.xlabel('amplicons (sorted)', size=14)
            plt.ylabel('# reads (log10)', size=14)
            plt.title(lib, size=14)
            plt.savefig('%s/%s_uniformity.png' % (uniformity_folder, lib), bbox_inches='tight')
            if to_plot:
                plt.show()
            else:
                plt.close()
        #------------------- uniformity -------------------#


        #------------------- dimer -------------------#
        if DIMER_SWITCH:
            print('analyzing dimers...')
            # load hash index for searching
            with open('%s/PrimerHeadDict_tube%d.json' % (index_folder, n_tube)) as f:
                primerhash = json.load(f)

            # count every unique dimer
            amp_cnt = Counter()
            with gzip.open('%s/%s_dimer.fastq.gz' % (dimer_folder, lib), 'rb') as f:
                for n, line in enumerate(f):
                    if n % 4 == 1:
                        # remove '\n' at the end
                        line = line.decode('ascii')[:-1]
                        amp_cnt[line] += 1

            # primer pair analysis
            df_amp, mtx_count, df_pair = primer_pairwise(primerhash, amp_cnt)
            N = (n+1)/4
            matched = sum(df_pair['count'])
            if verbose:
                print('Total dimer reads: %d' % N)
                print('matched reads: %d, %.3f%%' % (matched, 100*matched/N))
            df_amp.to_csv('%s/%s_dimer_list.csv' % (dimer_folder, lib))
            df_pair.to_csv('%s/%s_dimer_pairs.csv' % (dimer_folder, lib))

            # plot and save heatmap
            font_size = 20
            linewidth = 0.5
            mtx_count_log10 = np.log10(mtx_count+1)

            plt.figure(figsize=[12,10])
            sns.heatmap(mtx_count_log10, linewidth=linewidth)
            plt.title('%s dimer: %d, %.3f%%' % (lib, matched, 100*matched/N), fontsize=font_size)
            plt.xlabel('rP', fontsize=font_size)
            plt.ylabel('fP', fontsize=font_size)
            plt.savefig('%s/%s_dimer_pairs_log10.svg' % (dimer_folder, lib), bbox_inches='tight')
            if to_plot:
                plt.show()
            else:
                plt.close()
        #------------------- dimer -------------------#


        #------------------- non-specific -------------------#
        if NS_SWITCH:
            print('analyzing non-specific...')
            # load hash index for searching
            with open('%s/PrimerHeadDict_tube%d.json' % (index_folder, n_tube)) as f:
                primerhash = json.load(f)

            # count each unique non-specific amplicon
            amp_cnt = Counter()
            N = 0
            samfile = pysam.AlignmentFile(bam_path, 'rb')
            for record in samfile.fetch(until_eof=True):
                if not record.is_unmapped:
                    continue
                N += 1
                seq = record.query_sequence # read sequence
                amp_cnt[seq] += 1

            # primer pair analysis
            df_amp, mtx_count, df_pair = primer_pairwise(primerhash, amp_cnt)
            matched = sum(df_pair['count'])
            if verbose:
                print('Total non-specific reads: %d' % N)
                print('matched reads: %d, %.3f%%' % (matched, 100*matched/N))
            df_amp.to_csv('%s/%s_NS_list.csv' % (NS_folder, lib))
            df_pair.to_csv('%s/%s_NS_pairs.csv' % (NS_folder, lib))

            # plot and save heatmap
            font_size = 20
            linewidth = 0.5
            mtx_count_log10 = np.log10(mtx_count+1)
            
            plt.figure(figsize=[12,10])
            sns.heatmap(mtx_count_log10, linewidth=linewidth)
            plt.title('%s non-specific: %d, %.3f%%' % (lib, matched, 100*matched/N), fontsize=font_size)
            plt.xlabel('rP', fontsize=font_size)
            plt.ylabel('fP', fontsize=font_size)
            plt.savefig('%s/%s_NS_pairs_log10.svg' % (NS_folder, lib), bbox_inches='tight')
            if to_plot:
                plt.show()
            else:
                plt.close()
        #------------------- non-specific -------------------#


        #------------------- overall metrics -------------------#
        if TRIM_ALIGN_SWITCH and UNIFORMITY_SWITCH:
            all_total_arr = np.array(all_total)
            all_noerror_arr = np.array(all_noerror)
            all_short_arr = np.array(all_short)
            all_aligned_arr = np.array(all_aligned)

            pass_th = all_noerror_arr-all_short_arr

            noerror_rate = all_noerror_arr/all_total_arr
            pass_th_rate = pass_th/all_noerror_arr
            trim_rate = pass_th/all_total_arr
            align_rate = all_aligned_arr/pass_th
            on_target_rate = all_aligned_arr/all_total_arr

            df = pd.DataFrame({
                'lib_name': lib_list[:i_lib+1], 
                'total_reads': all_total_arr, 
                'unmerged': all_unmerged, 
                'unmatched': all_unmatched, 
                'search_succeeded': all_noerror_arr, 
                'search_succeeded_rate': noerror_rate, 
                'pass_lenth_th': pass_th, 
                'pass_lenth_th_rate': pass_th_rate, 
                'trim_rate': trim_rate, 
                'aligned': all_aligned_arr, 
                'align_rate': align_rate, 
                'on-target_rate': on_target_rate, 
                'MAD': all_mad, 
                '(max-min)/max': all_max_min, 
                '(95-5)/95': all_95_5, 
                'n_dropout': all_n_dropout, 
                'med_UMI_count': all_med_UMI
            })
            df.to_excel('metrics.xlsx')
        #------------------- overall metrics -------------------#


        # save reads per plex
        if UNIFORMITY_SWITCH:
            with pd.ExcelWriter('reads_per_amplicon.xlsx') as writer:
                for i in range(total_tube):
                    df = pd.DataFrame(read_dict[i])
                    df.to_excel(writer, sheet_name='tube_%d' % (i+1), index=True)

        # save UMI number
        if UMI_SWITCH:
            with pd.ExcelWriter('UMI_number.xlsx') as writer:
                for i in range(total_tube):
                    df = pd.DataFrame(umi_num_dict[i])
                    df.to_excel(writer, sheet_name='tube_%d' % (i+1), index=True)
            with pd.ExcelWriter('downsampling_CV.xlsx') as writer:
                for i in range(total_tube):
                    df = pd.DataFrame(DS_cv_dict[i])
                    df.to_excel(writer, sheet_name='tube_%d' % (i+1), index=True)


        if not keep_sam:
            # delete sam files
            os.system('rm %s' % sam_path)


        # end of loop
        print('%s finished in %.2fs\n' % (lib, time()-tik))

    print('finished!')


def main(config):
    """
    Example of config: 
    config = {
        # e.g., 1_ABC-pippin_S3 correspond to 
        # 1_ABC-pippin_S3_L001_R1_001.fastq.gz + 1_ABC-pippin_S3_L001_R2_001.fastq.gz
        'lib_list': [
            '1-SETD2NSD2-0p2-SNV-CW-5_S5', 
            '1-SETD2NSD2-0p5-SNV-KN-9_S9', 
            '2-SETD2NSD2-0p5-SNV-KN-10_S10'
        ], 
        
        # panel information
        'design_file_path': '../EpizymeLongAmp_design_redesign20210122.xlsx', # path to the design spread sheet
        'panel_type': 'long', # 'long': long amplicon panel; 'short': short amplicon panel
        'total_tube': 2, # number of tubes
        
        # fastq path format (NOTE: need to change format if merged from multiple lanes)
        'r1_path_format': 'raw_reads/%s_R1.fastq.gz',  # %s is library name
        'r2_path_format': 'raw_reads/%s_R2.fastq.gz', 
        
        # output
        'verbose': True, # print in notebook
        'to_plot': False, # show plots in notebook
        'keep_sam': False, # keep sam files
        'save_raw_UMI': False, # save UMI grouping file
        
        # performace
        'performace_level': 2, # 1 to use fewer RAM and half of all logic cores, 2 to use more RAM and all cores
        # e.g., set to 1 when running on laptop or 2 when running on AWS or desktop
        
        # function block switches
        'TRIM_ALIGN_SWITCH': True, # adapter trim & bowtie2 alignment (necessary for other function blocks)
        'UMI_SWITCH': True, # UMI grouping, UMI voting and variant call (need bam files)
        'UNIFORMITY_SWITCH': True, # uniformity analysis and plots (need bam files)
        'DIMER_SWITCH': False, # dimer analysis and plots (need dimer.fastq.gz)
        'NS_SWITCH': False # non-specific analysis and plots (need bam files)
    }
    """
    with open('config.json', 'w') as f:
        json.dump(config, f, indent=4)

    if config['BUILD_INDEX']:
        build_index(config['design_file_path'], config['total_tube'])
    del config['design_file_path'], config['BUILD_INDEX']
    run_pipeline(**config)
