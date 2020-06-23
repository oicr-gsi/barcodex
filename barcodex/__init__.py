# -*- coding: utf-8 -*-
"""
Created on Tue May 12 14:23:08 2020

@author: rjovelin
"""

import gzip
import os
from itertools import zip_longest
import regex
import argparse
import json
import time

def _is_gzipped(filename):
    '''
    (str) -> bool

    Returns True if the file is gzipped and False otherwise

    Parameters
    ----------
    
     - filename (str): File name or file path    
    '''
    
    # open file in rb mode
    infile = open(filename, 'rb')
    header = infile.readline()
    infile.close()
    if header.startswith(b'\x1f\x8b\x08'):
        return True
    else:
        return False


def _open_fastq(fastq):
    '''
    (str) -> _io.TextIOWrapper
    
    Returns an open fastq file with file handler at the begining of the file
    
    Parameters
    ----------
    
    - fastq (str): Path to the fastq file (compressed or not)
    '''
    
    # open input fastq
    if _is_gzipped(fastq) == True:
        infile = gzip.open(fastq, 'rt')
    else:
        infile = open(fastq)
    return infile


def _add_umi_to_readname(readname, UMI, separator):
    '''
    (str, str, str) -> str
    
    Returns the read name with the UMI sequence separated by separator
    
    Parameters
    ----------
    
    - readname (str): Read header
    - UMI (str): UMI sequence
    - separator (str): String separating the UMI sequence and part of the read header

    Examples
    --------
    >>> _add_umi_to_readname('@MISEQ753:114:000000000-D6365:1:1101:12254:19531 1:N:0:ATCACG', 'ATCG', '_')
    '@MISEQ753:114:000000000-D6365:1:1101:12254:19531_ATCG 1:N:0:ATCACG'
    >>> _add_umi_to_readname('@MISEQ753:114:000000000-D6365:1:1101:12254:19531 1:N:0:ATCACG', 'ATCGAT', ';')
    '@MISEQ753:114:000000000-D6365:1:1101:12254:19531;ATCGAT 1:N:0:ATCACG'
    '''
    
    readname = readname.split(' ')
    readname[0] = readname[0] + separator + UMI
    readname = ' '.join(readname)
    return readname


def _is_pattern_sequence(pattern):
    '''
    (str) -> bool
    
    Returns True if all elements of pattern are valid nucleotides
    
    parameters
    ----------
    - pattern (str): Pattern to be extracted from reads
    
    Examples
    --------
    >>> _is_pattern_sequence('(?<umi_1>.{3})AA')
    False
    >>> _is_pattern_sequence('ATCG')
    True
    >>> _is_pattern_sequence('ATCGNNNAX')
    False
    '''
    
    return all(map(lambda x: x in 'atcgnATCGN', set(pattern)))
    

def _find_pattern_umi(pattern):
    '''
    (str) -> list

    Returns a list with UMI sequence and eventually the spacer sequence

    Parameters
    ----------
    - pattern (str): String sequence used for matching and extracting UMis from reads.
                     Must look like NNNATCG or NNN. UMI nucleotides are labeled with "N".
                     Spacer nucleotides following Ns are used for matching UMIs but are
                     discarded from reads    
    Examples
    --------
    >>> _find_pattern_umi('NNNNATCG')
    ['NNNN', 'ATCG']
    >>> _find_pattern_umi('NNNNATCGNNN')
    ['NNNN', 'ATCG', 'NNN']
    >>> _find_pattern_umi('ATCGNNN')
    ['ATCG', 'NNN']
    >>> _find_pattern_umi('ATCG')
    ['ATCG']
    >>> _find_pattern_umi('NNNN')
    ['NNNN']
    '''

    # separate UMI and spacer from pattern
    if (len(set(pattern)) == 1 and list(set(pattern))[0] == 'N') or pattern == '':
        L = [pattern]
    else:
        L = list(map(lambda x: 'N' if x == '' else x, pattern.split('N')))
    # initiate list and seq 
    P, s = [], ''
    for i in range(len(L)):
        if L[i] == 'N':
            s += L[i]
        else:
            if s != '':
                P.append(s)
            P.append(L[i])
            s = ''
        if i == len(L)-1:
            if s != '':
                P.append(s)
    return P        
    

def _check_pattern_sequence(pattern):
    '''
    (str) -> None
   
    Raise a ValueError if the string pattern does not look like NNN or NNNATCG
    
    Parameters
    ----------
    - pattern (str): String sequence used for matching and extracting UMis from reads.
                     Must look like NNNATCG or NNN. UMI nucleotides are labeled with "N".
                     Spacer nucleotides following Ns are used for matching UMIs but are
                     discarded from reads    
    
    Examples
    --------
    
    >>> _check_pattern_sequence('NNNatcgatc')
    >>> _check_pattern_sequence('NNnnatcgatc')
    ValueError: String pattern must look like NNNNATCG or NNNN
    >>> _check_pattern_sequence('NNNNatcATCG')
    >>> _check_pattern_sequence('NNNNatcATCGNNNN')
    ValueError: String pattern must look like NNNNATCG or NNNN
    >>> _check_pattern_sequence('atcATCGNNNN')
    ValueError: String pattern must look like NNNNATCG or NNNN
    '''
    
    P = _find_pattern_umi(pattern)
    if len(P) > 2 or len(P) == 0:
        raise ValueError('String pattern must look like NNNNATCG or NNNN')
    else:
        if not(len(set(P[0])) == 1 and list(set(P[0]))[0] == 'N'):
            raise ValueError('String pattern must look like NNNNATCG or NNNN')
        if len(P) == 2:
            if not all(map(lambda x: x in 'atcgATCG', set(P[1]))):
                raise ValueError('String pattern must look like NNNNATCG or NNNN')


def _check_extraction_mode(pattern1, pattern2):
    '''
    (str | None, str | None) -> None
    
    Raise a ValueError if pattern and pattern2 are not both a string sequence or a regex
    
    Parameters
    ----------
    - pattern1 (str or None): String sequence or regular expression used for matching and extracting UMis from reads in FASTQ 1.
                             None if UMIs are extracted only from FASTQ 2 
    - pattern2 (str or None): String sequence or regular expression used for matching and extracting UMis from reads in FASTQ 2.
                              None if UMIs are extracted only from FASTQ 1.
                              
    Examples
    --------    
    >>> _check_extraction_mode('NNNATCG', 'NNNNGTCG')
    >>>  _check_extraction_mode('NNNATCG', '(?<umi_1>.{3})AA')
    ValueError: Both patterns must be either string sequences or regex
    >>> _check_extraction_mode('(?<discard_1>.+)(?<umi_1>.{3})(?discard_2>TT)', '(?<umi_1>.{3})AA')
    '''
    
    if pattern1 and pattern2:
        if _is_pattern_sequence(pattern1) == True and _is_pattern_sequence(pattern2) == False:
            raise ValueError('Both patterns must be either string sequences or regex')
        elif _is_pattern_sequence(pattern1) == False and _is_pattern_sequence(pattern2) == True:
            raise ValueError('Both patterns must be either string sequences or regex')


def _get_umi_spacer(pattern):
    '''
    (str) -> (str, str)
    
    Returns the UMI and spacer sequences in pattern
    
    Parameters
    ----------
    - pattern (str): String sequence used for matching and extracting UMis from reads.
                     Must look like NNNATCG or NNN. UMI nucleotides are labeled with "N".
                     Spacer nucleotides following Ns are used for matching UMIs but are
                     discarded from reads    
    
    Examples
    --------
    >>> _get_umi_spacer('NNNNAAAA')
    ('NNNN', 'AAAA')
    >>> _get_umi_spacer('NNNN')
    ('NNNN', '')
    >>> _get_umi_spacer('AAATCGC')
    ('AAATCGC', '')
    >>> _get_umi_spacer('AAATCGCNNN')
    ('AAATCGC', 'NNN')
    >>> _get_umi_spacer('NNNAAATCGCNNN')
    ('NNN', 'AAATCGC')
    '''
        
    P = _find_pattern_umi(pattern)
    if len(P) == 1:
        UMI, spacer = P[0], ''    
    else:
        UMI, spacer = P[0], P[1]
    return UMI, spacer
    

def _extract_from_sequence(read, UMI, spacer):
    '''
    (list, str, str) -> (str, str, str, str, str)
    
    Returns a tuple with the read sequence and qualities after barcode extraction,
    the umi sequence, the read sequence and qualities extracted from read.
    Or a tuple with empty strings when there is no match
    
    Parameters
    ----------
    - read (list): List of 4 strings from a single read
    - UMI (str): UMI nucleotides are labeled with "N" (eg, NNNN)
    - spacer (str): Spacer sequence following the UMI. Can be the empty string or
                    any nucleotides from 'ATCGatcg'. Spacer sequences are extracted
                    and discarded from reads 
    
    Examples
    --------
    read = ['@MISEQ753:39:000000000-BDH2V:1:1101:17521:1593 1:N:0:', 'TCATGTCTGCTAATGGGAAAGAGTGTCCTAACTGTCCCAGATCGTTTTTTCTCACGTCTTTTCTCCTTTCACTTCTCTTTTTCTTTTTCTTTCTTCTTCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT', '+',  '1>1A1DDF11DBDGFFA111111D1FEEG31AD1DAA1110BA00000//01A2A/B/B/212D2111D1222D12122B1B01D1@101112@D2D12BB##################################################']
    >>> _extract_from_sequence(read, 'NNNNNNNNNNNN', 'ATGGGAAAGAGTGTCC')
    ('TAACTGTCCCAGATCGTTTTTTCTCACGTCTTTTCTCCTTTCACTTCTCTTTTTCTTTTTCTTTCTTCTTCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT',
    'G31AD1DAA1110BA00000//01A2A/B/B/212D2111D1222D12122B1B01D1@101112@D2D12BB##################################################',
    'TCATGTCTGCTA',
    'TCATGTCTGCTAATGGGAAAGAGTGTCC',
    '1>1A1DDF11DBDGFFA111111D1FEE')
    >>> _extract_from_sequence(read, 'NNNNNNNNNN', 'ATGGCATCG')
    ('', '', '', '', '')
    >>> _extract_from_sequence(read, 'NNNNNNNNNN', '')
    ('TAATGGGAAAGAGTGTCCTAACTGTCCCAGATCGTTTTTTCTCACGTCTTTTCTCCTTTCACTTCTCTTTTTCTTTTTCTTTCTTCTTCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT',
     'DBDGFFA111111D1FEEG31AD1DAA1110BA00000//01A2A/B/B/212D2111D1222D12122B1B01D1@101112@D2D12BB##################################################',
     'TCATGTCTGC',
     'TCATGTCTGC',
     '1>1A1DDF11')
    '''
    
    # initialize variables
    seq, qual, umi_seq, extracted_seq, extracted_qual = '', '', '', '', ''
    
    # extract UMI starting at begining of the read
    if spacer in read[1]:
        if spacer == read[1][len(UMI): len(UMI) + len(spacer)]:
            umi_seq = read[1][: len(UMI)]
            extracted_seq = read[1][: len(UMI) + len(spacer)]
            extracted_qual = read[3][: len(UMI) + len(spacer)]
            seq = read[1][len(UMI) + len(spacer):]
            qual = read[3][len(UMI) + len(spacer):]
    
    return seq, qual, umi_seq, extracted_seq, extracted_qual
    

def _extract_from_regex(read, p, full_match=False):
    '''
    (list, _regex.Pattern, bool) -> (str, str, str, str, str)

    Returns a tuple with the read sequence and qualities after barcode extraction,
    the umi sequence, the read sequence and qualities extracted from read
    Or a tuple with empty strings when there is no match
    
    Parameters
    ----------
    - read (list): List of 4 strings from a single read
    - p (_regex.Pattern): Compiled regex pattern used for matching pattern in read sequence
    - full_match (bool): True if the regular expression needs to match the entire read sequence 
    
    Examples
    --------
    read = ['@MISEQ753:39:000000000-BDH2V:1:1101:17521:1593 1:N:0:', 'TCATGTCTGCTAATGGGAAAGAGTGTCCTAACTGTCCCAGATCGTTTTTTCTCACGTCTTTTCTCCTTTCACTTCTCTTTTTCTTTTTCTTTCTTCTTCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT', '+',  '1>1A1DDF11DBDGFFA111111D1FEEG31AD1DAA1110BA00000//01A2A/B/B/212D2111D1222D12122B1B01D1@101112@D2D12BB##################################################']
    # UMI starts at the beginning of the read
    >>> _extract_from_regex(read, regex.compile('(?<umi_1>.{12})(?<discard_1>ATGGGAAAGAGTGTCC)'), full_match=False)
    ('TAACTGTCCCAGATCGTTTTTTCTCACGTCTTTTCTCCTTTCACTTCTCTTTTTCTTTTTCTTTCTTCTTCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT',
     'G31AD1DAA1110BA00000//01A2A/B/B/212D2111D1222D12122B1B01D1@101112@D2D12BB##################################################',
     'TCATGTCTGCTA',
     'TCATGTCTGCTAATGGGAAAGAGTGTCC',
     '1>1A1DDF11DBDGFFA111111D1FEE')
    # match the entire read sequence
    >>> _extract_from_regex(read, regex.compile('(?<umi_1>.{12})(?<discard_1>ATGGGAAAGAGTGTCC)'), full_match=True)
    ('', '', '', '', '')
    # contruct the regex to match the entire read sequence
    >>> _extract_from_regex(read, regex.compile('(?<umi_1>.{12})(?<discard_1>ATGGGAAAGAGTGTCC)[ATCG]*'), full_match=True)
    ('TAACTGTCCCAGATCGTTTTTTCTCACGTCTTTTCTCCTTTCACTTCTCTTTTTCTTTTTCTTTCTTCTTCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT',
     'G31AD1DAA1110BA00000//01A2A/B/B/212D2111D1222D12122B1B01D1@101112@D2D12BB##################################################',
     'TCATGTCTGCTA',
     'TCATGTCTGCTAATGGGAAAGAGTGTCC',
     '1>1A1DDF11DBDGFFA111111D1FEE')
    
    # UMI does not start at the beginning of the read
    read = ['@MISEQ753:39:000000000-BDH2V:1:1101:17521:1593 1:N:0:', 'ATCATGTCTGCTAATGGGAAAGAGTGTCCTAACTGTCCCAGATCGTTTTTTCTCACGTCTTTTCTCCTTTCACTTCTCTTTTTCTTTTTCTTTCTTCTTCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT', '+', 'B1>1A1DDF11DBDGFFA111111D1FEEG31AD1DAA1110BA00000//01A2A/B/B/212D2111D1222D12122B1B01D1@101112@D2D12BB##################################################']    
    # first nucleotide is part of the new read sequence
    >>> _extract_from_regex(read, regex.compile('(?<umi_1>.{12})(?<discard_1>ATGGGAAAGAGTGTCC)'), full_match=False)
    ('ATAACTGTCCCAGATCGTTTTTTCTCACGTCTTTTCTCCTTTCACTTCTCTTTTTCTTTTTCTTTCTTCTTCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT',
     'BG31AD1DAA1110BA00000//01A2A/B/B/212D2111D1222D12122B1B01D1@101112@D2D12BB##################################################',
     'TCATGTCTGCTA',
     'TCATGTCTGCTAATGGGAAAGAGTGTCC',
     '1>1A1DDF11DBDGFFA111111D1FEE')
    # force UMI to start at the beginning of the read sequence
    >>> _extract_from_regex(read, regex.compile('(?<umi_1>^.{12})(?<discard_1>ATGGGAAAGAGTGTCC)'), full_match=False)
    ('', '', '', '', '')
    # discard nuceotides upstream of UMI
    >>> _extract_from_regex(read, regex.compile('(?<discard_1>.*)(?<umi_1>.{12})(?<discard_2>ATGGGAAAGAGTGTCC)'), full_match=False)
    ('TAACTGTCCCAGATCGTTTTTTCTCACGTCTTTTCTCCTTTCACTTCTCTTTTTCTTTTTCTTTCTTCTTCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT',
     'G31AD1DAA1110BA00000//01A2A/B/B/212D2111D1222D12122B1B01D1@101112@D2D12BB##################################################',
     'TCATGTCTGCTA',
     'ATCATGTCTGCTAATGGGAAAGAGTGTCC',
     'B1>1A1DDF11DBDGFFA111111D1FEE')
    '''
    
    # initialize variables
    seq, qual, umi_seq, extracted_seq, extracted_qual = '', '', '', '', ''
        
    # look for a match in read sequence
    if full_match == False:
        # scan through the string looking for a match
        m = p.search(read[1])
    elif full_match == True:
        # match if the whole string matches pattern 
        m = p.fullmatch(read[1])
    # process if match is found
    if m:
        # collect umi, discard positions
        umi_pos, discard_pos = [], []
        for i in m.groupdict():
            if 'umi' in i:
                umi_pos.append(m.span(i))
            elif 'discard' in i:
                discard_pos.append(m.span(i))
        # sort umi and discard positions
        umi_pos.sort()
        discard_pos.sort()
        # get umi sequences
        umi_seq = ''.join([read[1][i[0]:i[1]] for i in umi_pos])
        # get indices of extracted sequences
        extracted_pos = sorted(umi_pos + discard_pos)
        # get indices of remaining sequence
        removed = [i for j in umi_pos + discard_pos for i in list(range(j[0], j[1]))]
        remaining = sorted([i for i in range(len(read[1])) if i not in removed])
        
        # get extracted sequence and qualities
        extracted_seq = ''.join([read[1][i[0]: i[1]] for i in extracted_pos])
        extracted_qual = ''.join([read[3][i[0]: i[1]] for i in extracted_pos])
        
        # get read seq and qual after extraction
        seq = ''.join([read[1][i] for i in remaining])
        qual = ''.join([read[3][i] for i in remaining])
           
    return seq, qual, umi_seq, extracted_seq, extracted_qual


def _get_read(fastq_file):
    """
    (_io.TextIOWrapper) -- > itertools.zip_longest
   
    Returns an iterator slicing the fastq into 4-line reads.
    Each element of the iterator is a tuple containing read information
    
    Parameters
    ----------
    
    - fastq_file (_io.TextIOWrapper): Fastq file opened for reading in plain text mode
    """
    args = [iter(fastq_file)] * 4
    return zip_longest(*args, fillvalue=None)


def _read_cleanup(read):
    '''
    (tuple) -> list
    
    Takes a tuple with tuple(s) containing read information and retutn a list
    of lists in which all elements are stripped of trailing characters
    
    Parameters
    ----------
    - read (tuple): A tuple with one or more tuples, each containing 4 strings from a read    
    
    Examples
    -------
    read = (('@M00146:137:000000000-D7KWF:1:1102:19596:10317 1:N:0:GTTCTCGT\n',
             'ACTGTTGAGATACTTAGTAATAAATTAAATAAACATTTCTAAAAGAGTATTCTACATTTTTAGCCTAAACATATAAGAGAAAGCATCTGAAGCAGTCATGTCACACAGTAGAGATAATTGTTGATGATGAAATAATCACAGTAGAGGTCAT\n',
             '+\n',
             'CCCCCFFFFCFFGGGGGGGGGGHHHHHHHHHHHHHHHHHHHGGHHGHGHHHHHHHHHIHHHGHIHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHGHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH'),
            ('@M00146:137:000000000-D7KWF:1:1102:19596:10317 2:N:0:GTTCTCGT\n',
             'CCCATGACCTCTACTGTGATTATTTCATCATCAACAATTATCTCTACTGTGTGACATGACTGCTTCAGATGCTTTCTCTTATATGTTTAGGCTAAAAATGTAGAATACTCTTTTAGAAATGTTTATTTAATTTATTACTAAGTATCTCAAC\n',
             '+\n',
             'BCBCCFFFFFFFGGGGGGGGGGHHHHHHHHHHHHHHGHHHHHHHHHHHHGHHHHHFHHHHHHHHHHHHGIHHHHHHHHHHGHHHHHHHHGHGHHHHHHHHHHHHHHHHHHGHHHHHGHGHHHHHHHHHHHHHHHHHIHHHHHHHHHHHHHH'))
    
    >>> _read_cleanup(read)
    [['@M00146:137:000000000-D7KWF:1:1102:19596:10317 1:N:0:GTTCTCGT',
      'ACTGTTGAGATACTTAGTAATAAATTAAATAAACATTTCTAAAAGAGTATTCTACATTTTTAGCCTAAACATATAAGAGAAAGCATCTGAAGCAGTCATGTCACACAGTAGAGATAATTGTTGATGATGAAATAATCACAGTAGAGGTCAT',
      '+',
      'CCCCCFFFFCFFGGGGGGGGGGHHHHHHHHHHHHHHHHHHHGGHHGHGHHHHHHHHHIHHHGHIHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHGHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH'],
     ['@M00146:137:000000000-D7KWF:1:1102:19596:10317 2:N:0:GTTCTCGT',
      'CCCATGACCTCTACTGTGATTATTTCATCATCAACAATTATCTCTACTGTGTGACATGACTGCTTCAGATGCTTTCTCTTATATGTTTAGGCTAAAAATGTAGAATACTCTTTTAGAAATGTTTATTTAATTTATTACTAAGTATCTCAAC',
      '+',
      'BCBCCFFFFFFFGGGGGGGGGGHHHHHHHHHHHHHHGHHHHHHHHHHHHGHHHHHFHHHHHHHHHHHHGIHHHHHHHHHHGHHHHHHHHGHGHHHHHHHHHHHHHHHHHHGHHHHHGHGHHHHHHHHHHHHHHHHHIHHHHHHHHHHHHHH']]
    '''
    
    read = list(map(lambda x: list(x), read))
    for i in range(len(read)):
        read[i] = list(map(lambda x: x.strip(), read[i]))
    return read


def _check_fastq_sync(L):
    '''
    (list) -> None
    
    Raise a ValueError if the read headers in L are not from the same paired reads
    (ie, if the fastqs the reads originate from are not synced)
    
    Parameters
    ----------
    
    - a list of read headers
    
    Examples
    --------
    >>> _check_fastq_sync(['@MISEQ753:114:000000000-D6365:1:1101:15443:1350 1:N:0:ATCACG\n',
    '@MISEQ753:114:000000000-D6365:1:1101:15443:1350 2:N:0:ATCACG\n'])
    >>> _check_fastq_sync(['@MISEQ753:114:000000000-D6378:1:1102:15450:1350 1:N:0:ATCGGG\n',
    '@MISEQ753:114:000000000-D6365:1:1101:15443:1350 2:N:0:ATCACG\n'])
    ValueError: Fastqs are not synced
    '''
    
    readnames = []
    for i in L:
        readnames.append(i.split(' ')[0])
    if len(set(readnames)) > 1:
        raise ValueError('Fastqs are not synced')
        
    
def _check_input_output(r1_in, r1_out, data='single', inline_umi=True,
                        r2_in=None, r2_out=None, r3_in=None):
    '''
    (list, str, str, bool, list | None, str | None, list | None) -> None
    
    Raises a ValueError if input / output fastqs are not compatible with paired
    or single end sequencing data and with inline UMIs or UMIs in separate fastq.
    
    Parameters
    ----------
    
    - r1_in (list): Path(s) to the input FASTQ 1 (compressed or not) 
    - r1_out (str): Path to the output FASTQ 1 with reads re-headered with UMI sequence
    - data (str): Indicates if single or paired end sequencing data
    - inline_umi (bool): True if UMIs are inline with reads and False otherwise
    - r2_in (list or None): Path(s) to the input FASTQ 2 (compressed or not)
    - r2_out (str or None): Path to the output FASTQ 2 with reads re-headered with UMI sequence
    - r3_in (list or None): Path(s) to input FASTQ 3 for paired end sequences with non-inline UMIs
    
    Examples
    --------
    >>> _check_input_output(['infile_1.fastq'], 'outputfile_1.fastq', data='single', inline_umi=True, r2_in=None, r2_out=None, r3_in=None)
    >>> _check_input_output(['infile_1.fastq', 'infile_2.fastq'], 'outputfile_1.fastq', data='single', inline_umi=True, r2_in=None, r2_out=None, r3_in=None)
    >>> _check_input_output(['infile_1.fastq'], 'outputfile_1.fastq', data='single', inline_umi=True, r2_in=None, r2_out=None, r3_in=['infile_3.fastq'])
    ValueError: Expecting single end sequences with inline UMIs. Paths to r1 I/O fastqs required. Paths to r2 I/O fastqs and r3 input fastq not needed
    >>> _check_input_output(['infile_1.fastq'], 'outputfile_1.fastq', data='single', inline_umi=False, r2_in=None, r2_out=None, r3_in=['infile_3.fastq'])
    ValueError: Expecting single end sequences with out of read UMIs. Paths to r1 I/O and r2 input fastqs required
    >>> _check_input_output(['infile_1.fastq'], 'outputfile_1.fastq', data='single', inline_umi=False, r2_in=['infile_2.fastq'], r2_out=None, r3_in=['infile_3.fastq'])
    ValueError: Expecting single end sequences with out of read UMIs. Paths to r2 output and r3 input fastq not needed
    >>> _check_input_output(['infile_1.fastq'], 'outputfile_1.fastq', data='single', inline_umi=False, r2_in=['infile_2.fastq'], r2_out=None, r3_in=None)
    >>> _check_input_output(['infile_1.fastq'], 'outputfile_1.fastq', data='single', inline_umi=True, r2_in=['infile_2.fastq'], r2_out=None, r3_in=None)
    ValueError: Expecting single end sequences with inline UMIs. Paths to r1 I/O fastqs required. Paths to r2 I/O fastqs and r3 input fastq not needed
    >>> _check_input_output(['infile_1.fastq'], 'outputfile_1.fastq', data='paired', inline_umi=True, r2_in=['infile_2.fastq'], r2_out='outputfile_2.fastq', r3_in=None)
    >>> _check_input_output('infile_1.fastq', 'outputfile_1.fastq', data='paired', inline_umi=True, r2_in='infile_2.fastq', r2_out='outputfile_2.fastq', r3_in='inputfile_3.fastq')
    ValueError: Expecting paired end sequences with inline UMIs. Paths to r1 and r2 I/O fastqs required. Path to r3 input not needed
    '''
    
    if data == 'paired':
        if inline_umi:
            # requires r1 and r2 but not r3
            # requires pattern, pattern2 optional
            if any(map(lambda x: x is None, [r1_in, r1_out, r2_in, r2_out])):
                raise ValueError('Expecting paired end sequences with inline UMIs. Paths to r1 and r2 I/O fastqs required')
            if not r3_in is None:
                raise ValueError('Expecting paired end sequences with inline UMIs. Paths to r1 and r2 I/O fastqs required. Path to r3 input not needed')
        elif not inline_umi:
            # requires r1, r2 and r3
            # requires pattern1, pattern2 not needed
            if any(map(lambda x: x is None, [r1_in, r1_out, r2_in, r2_out, r3_in])):
                raise ValueError('Expecting paired end sequences with out of read UMIs. Paths to r1 and r2 I/O fastqs and to r3 input fastq required')
    elif data == 'single':
        if inline_umi:
            # requires r1, but not r2 and not r3
            if any(map(lambda x: x is None, [r1_in, r1_out])):
                raise ValueError('Expecting single end sequences with inline UMIs. Paths to r1 I/O fastqs required')
            if any(map(lambda x: x is not None, [r2_in, r2_out, r3_in])):
                raise ValueError('Expecting single end sequences with inline UMIs. Paths to r1 I/O fastqs required. Paths to r2 I/O fastqs and r3 input fastq not needed')
        elif not inline_umi:
            # requires r1 and r2, r3 not needed
            if any(map(lambda x: x is None, [r1_in, r1_out, r2_in])):
                raise ValueError('Expecting single end sequences with out of read UMIs. Paths to r1 I/O and r2 input fastqs required')
            if any(map(lambda x: x is not None, [r2_out, r3_in])):
                raise ValueError('Expecting single end sequences with out of read UMIs. Paths to r2 output and r3 input fastq not needed')


def _check_input_files(L1, L2, L3):
    '''
    (list, list | None, list | None)

    Raises a ValueError if different number of input files for paired data
        
    Parameters
    ----------
    - L1 (list): List of read1 input files
    - L2 (list or None): List of read2 input files if paired data or None
    - L3 (list or None): List of read3 input files for paired data with non-inline UMIs or None

    Examples
    --------
    >>> _check_input_files(['r1.1', 'r1.2'], ['r2.1', 'r2.2'], [])
    >>> _check_input_files(['r1.1', 'r1.2'], ['r2.1', 'r2.2'], ['r3.1', 'r3.2'])
    >>> _check_input_files(['r1.1', 'r1.2'], ['r2.1', 'r2.2'], ['r3.1'])
    ValueError: Expecting same number of input files
    >>> _check_input_files(['r1.1', 'r1.2'], ['r2.1'], [])
    ValueError: Expecting same number of input files   
    '''

    if L2:
        if len(L1) != len(L2):
            raise ValueError('Expecting same number of input files')
    if L3:
        if len(L1) != len(L3):
            raise ValueError('Expecting same number of input files')


def _group_input_files(L1, L2, L3):
    '''
    (list, list | None, list | None) -> list
    
    Returns a list of same size tuples, each containing the paths of read1, read2
    and read3 fastqs when defined or None
        
    Parameters
    ----------
    - L1 (list): List of read1 input files
    - L2 (list or None): List of read2 input files if paired data or None
    - L3 (list or None): List of read3 input files for paired data with non-inline UMIs or None

    Examples
    --------
    >>> _group_input_files(['r1.1', 'r1.2'], ['r2.1', 'r2.2'], ['r3.1', 'r3.2'])
    [('r1.1', 'r2.1', 'r3.1'), ('r1.2', 'r2.2', 'r3.2')]
    >>> _group_input_files(['r1.1', 'r1.2'], ['r2.1', 'r2.2'], None)
    [('r1.1', 'r2.1', None), ('r1.2', 'r2.2', None)]
    >>> _group_input_files(['r1.1', 'r1.2'], None, None)
    [('r1.1', None, None), ('r1.2', None, None)]
    '''
    
    if L2 is None:
        L2 = [None] * len(L1)
    if L3 is None:
        L3 = [None] * len(L1)
    
    return list(zip(L1, L2, L3))
    

def _check_pattern_options(pattern1, pattern2=None, data='single', inline_umi=True):
    '''
    (str | None, str | None, str, bool) -> None
    
    Raise ValueError if pattern options are incompatible with single or paired
    sequencing data
    
    Parameters
    ----------
    
    - pattern1 (str or None): String sequence or regular expression used for matching and extracting UMis from reads in FASTQ 1.
                            None if UMIs are extracted only from FASTQ 2 
    - pattern2 (str or None): String sequence or regular expression used for matching and extracting UMis from reads in FASTQ 2.
                              None if UMIs are extracted only from FASTQ 1.
    - data (str): Indicates if single or paired end sequencing data
    - inline_umi (bool): True if UMIs are inline with reads and False otherwise
    
    Examples
    --------    
    >>> _check_pattern_options('NNNATCG', pattern2=None, data='single', inline_umi=False)
    >>> _check_pattern_options('NNNATCG', pattern2='ATCG', data='paired', inline_umi=False)
    ValueError: Expecting paired end sequences with UMIs in separate fastq. Requires pattern. Pattern2 is not needed
    >>> _check_pattern_options('NNNATCG', pattern2='ATCG', data='paired', inline_umi=True)
    >>> _check_pattern_options(None, pattern2='ATCG', data='paired', inline_umi=True)
    >>> _check_pattern_options(None, pattern2=None, data='paired', inline_umi=True)
    ValueError: Expecting paired end sequences with inline UMIs. At least 1 pattern is required
    >>> _check_pattern_options(None, pattern2='ATCG', data='paired', inline_umi=False)
    ValueError: Expecting paired end sequences with UMIs in separate fastq. Requires pattern. Pattern2 is not needed
    >>> _check_pattern_options(None, pattern2='ATCG', data='single', inline_umi=False)
    ValueError: Expecting single end sequences. Pattern required, pattern2 not needed
    '''
        
    if data == 'single':
        # expecting pattern but not pattern2
        if pattern2 or pattern1 is None:
            raise ValueError('Expecting single end sequences. Pattern1 required, pattern2 not needed')
    else:
        if inline_umi:
            # pattern1 is optional for paired end data with inline UMIs if pattern2 present
            if pattern1 is None and pattern2 is None:
                raise ValueError('Expecting paired end sequences with inline UMIs. At least 1 pattern is required')
        else:
            # pattern1 required for for non-inline UMI. pattern2 not needed
            if pattern1 is None or pattern2 is not None:
                raise ValueError('Expecting paired end sequences with UMIs in separate fastq. Requires pattern1. Pattern2 is not needed')
    
         
def _extract_umi_from_read(read, seq_extract, UMI, spacer, p, full_match):    
    '''
    (list, bool, str | None, str | None, _regex.Pattern | None, bool) -> (str, str, str, str, str)
    
    
    Returns a tuple with the read sequence and qualities after barcode extraction
    with a string pattern or a regex, the umi sequence, the read sequence and qualities
    extracted from read. Or a tuple with empty strings when there is no match
    
    Parameters
    ----------
    - read (list): List of 4 strings from a single read
    - UMI (str | None): UMI nucleotides are labeled with "N" (eg, NNNN)
    - spacer (str | None): Spacer sequence following the UMI. Can be the empty string or
                           any nucleotides from 'ATCGatcg'. Spacer sequences are extracted
                           and discarded from reads 
    - p (_regex.Pattern | None): Compiled regex pattern used for matching pattern in read sequence
    - full_match (bool): True if the regular expression needs to match the entire read sequence 
    
    Examples
    --------
    read = ['@MISEQ753:39:000000000-BDH2V:1:1101:17521:1593 1:N:0:', 'TCATGTCTGCTAATGGGAAAGAGTGTCCTAACTGTCCCAGATCGTTTTTTCTCACGTCTTTTCTCCTTTCACTTCTCTTTTTCTTTTTCTTTCTTCTTCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT', '+',  '1>1A1DDF11DBDGFFA111111D1FEEG31AD1DAA1110BA00000//01A2A/B/B/212D2111D1222D12122B1B01D1@101112@D2D12BB##################################################']
    >>> _extract_umi_from_read(read, True, 'NNNNNNNNNNNN', 'ATGGGAAAGAGTGTCC', None, True)
    ('TAACTGTCCCAGATCGTTTTTTCTCACGTCTTTTCTCCTTTCACTTCTCTTTTTCTTTTTCTTTCTTCTTCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT',
     'G31AD1DAA1110BA00000//01A2A/B/B/212D2111D1222D12122B1B01D1@101112@D2D12BB##################################################',
     'TCATGTCTGCTA',
     'TCATGTCTGCTAATGGGAAAGAGTGTCC',
     '1>1A1DDF11DBDGFFA111111D1FEE')
    >>> _extract_umi_from_read(read, True, 'NNNNNNNNNNNN', 'ATGGGAAAGAGTGTCC', regex.compile('(?P<umi_1>.{3})(?P<discard_1>.{2})'), True)
    ('TAACTGTCCCAGATCGTTTTTTCTCACGTCTTTTCTCCTTTCACTTCTCTTTTTCTTTTTCTTTCTTCTTCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT',
     'G31AD1DAA1110BA00000//01A2A/B/B/212D2111D1222D12122B1B01D1@101112@D2D12BB##################################################',
     'TCATGTCTGCTA',
     'TCATGTCTGCTAATGGGAAAGAGTGTCC',
     '1>1A1DDF11DBDGFFA111111D1FEE')
    >>> _extract_umi_from_read(read, False, 'NNNNNNNNNNNN', 'ATGGGAAAGAGTGTCC', regex.compile('(?<umi_1>.{12})(?<discard_1>ATGGGAAAGAGTGTCC)'), False)
    ('TAACTGTCCCAGATCGTTTTTTCTCACGTCTTTTCTCCTTTCACTTCTCTTTTTCTTTTTCTTTCTTCTTCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT',
     'G31AD1DAA1110BA00000//01A2A/B/B/212D2111D1222D12122B1B01D1@101112@D2D12BB##################################################',
     'TCATGTCTGCTA',
     'TCATGTCTGCTAATGGGAAAGAGTGTCC',
     '1>1A1DDF11DBDGFFA111111D1FEE')
    >>> _extract_umi_from_read(read, False, 'NNNNNNNNNNNN', 'ATGGGAAAGAGTGTCC', regex.compile('(?<umi_1>.{12})(?<discard_1>ATGGGAAAGAGTGTCC)'), True)
    ('', '', '', '', '')
    '''
    
    if seq_extract == True:
        # extraction using string sequence. assumes UMi starts at begining of read        
        L = _extract_from_sequence(read, UMI, spacer)
    else:
        L = _extract_from_regex(read, p, full_match)
    return L
    


def _get_read_patterns(pattern):
    '''
    (str) -> (bool, str | None, str | None, _regex.Pattern | None)
    
    Returns a tuple with pattern parameters for finding UMIs in read sequence:
    - a boolean indicating whether UMI is extracted with a string pattern or regex
    (or None if pattern is not defined)
    - the UMI sequence (labeled as Ns or None)
    - the spacer sequence (empty string, ATCGNatcgN nucleotides or None)
    - a regex pattern or None
    
    Parameters
    ----------
    
    - pattern (str): String sequence or regular expression used for matching and extracting UMis from reads

    Examples
    --------    
    >>> _get_read_patterns('NNNN')
    (True, 'NNNN', '', None)
    >>> _get_read_patterns('NNNNATCGA')
    (True, 'NNNN', 'ATCGA', None)
    >>> _get_read_patterns('NNNNATCGANNN')
    (True, 'NNNN', 'ATCGA', None)
    >>> _get_read_patterns('(?<umi_1>.{3})')
    (False, None, None, regex.Regex('(?<umi_1>.{3})', flags=regex.V0))
    >>> _get_read_patterns('(?<discard_1>.*)(?<umi_1>.{3})(?<discard_2>TT)')
    (False, None, None, regex.Regex('(?<discard_1>.*)(?<umi_1>.{3})(?<discard_2>TT)', flags=regex.V0))
    >>> _get_read_patterns(None)
    (None, None, None, None)
    >>> _get_read_patterns('')
    (None, None, None, None)
    '''
    
    # initialize variables
    seq_extract, UMI, spacer, p = None, None, None, None
    
    if pattern:
        # check if pattern is nucleotide string or regex
        if _is_pattern_sequence(pattern) == True:
            seq_extract = True
            # get UMi and spacer
            UMI, spacer = _get_umi_spacer(pattern)
        else:
            seq_extract = False
            # compile pattern
            p = regex.compile(pattern)

    return seq_extract, UMI, spacer, p


def _open_fastq_writing(output_file, compressed):
    '''
    (str, bool) -> _io.TextIOWrapper
    
    Returns a file handler for writing to output_file. The output_file is compressed
    with gzip (highest compression, level 9),  if compressed is True or is in 
    plain text file if False.
    
    Parameters
    ----------
    
    - output_file (str): Path to the output_file
    - compressed (bool): output_file is compressed with gzip if True 
    '''
    
    if compressed:
        newfile = gzip.open(output_file, 'wt')
    else:
        newfile = open(output_file, 'w')
    return newfile



def _get_files_extracted_reads(keep_extracted, data, inline_umi, pattern1, pattern2, r1_out, r2_out, compressed):
    '''
    (keep_extracted, data, inline_umi, pattern1, pattern2, r1_out, r2_out, compressed) -> (_io.TextIOWrapper | None, _io.TextIOWrapper | None, _io.TextIOWrapper | None)
    
    Returns a tuple with fastq files opened for writing extracted sequences 
    (UMIs and discarded sequences) if keep_discarded is True or a tuple with None if False.
    Each element of the tuple can be be a file handle for writing or None depending of arguments.
    Fastq files are written in the directory of r1_out and r2_out and are named
    by appending '.extracted_sequences.RN.fastq.gz' to r1_out and/or r2_out
    
    Parameters
    ----------
    
    - keep_extracted (bool): Write extracted sequences (UMIs and discarded sequences) to file if True
    - data (str): Indicates if single or paired end sequencing data
    - inline_umi (bool): True if UMIs are inline with reads and False otherwise
    - pattern1 (str or None): String sequence or regular expression used for matching and extracting UMis from reads in FASTQ 1
                             None if UMIs are extracted only from FASTQ 2 for paired end sequences
    - pattern2 (str or None): String sequence or regular expression used for matching and extracting UMis from reads in FASTQ 2
                              None if UMIs are extracted only from FASTQ 1 for paired end sequences
    - r1_out (str): Path to the output FASTQ 1 with reads re-headered with UMI sequence
    - r2_out (str or None): Path to the output FASTQ 2 with reads re-headered with UMI sequence    
                            None for single end read sequences
    - compressed (bool): output fastqs are compressed with gzip if True
    '''
    
    # initialize variables
    r1_extracted, r2_extracted, r3_extracted = None, None, None
    
    # add suffix to file name
    if compressed:
        suffix = '.extracted.{0}.fastq.gz'
    else:
        suffix = '.extracted.{0}.fastq'
    
    if keep_extracted:
        if inline_umi:
            if pattern1 is not None:
                r1_extracted = _open_fastq_writing(_remove_fastq_extension(r1_out) +  suffix.format('R1'), compressed)
            if pattern2 is not None:
                r2_extracted = _open_fastq_writing(_remove_fastq_extension(r2_out) +  suffix.format('R2'), compressed)
        else:
            if data == 'paired':
                r3_extracted = _open_fastq_writing(_remove_fastq_extension(r1_out) + suffix.format('R3'), compressed)
            elif data == 'single':
                r2_extracted = _open_fastq_writing(_remove_fastq_extension(r1_out) + suffix.format('R2'), compressed)

    return r1_extracted, r2_extracted, r3_extracted



def _get_files_discarded_reads(data, keep_discarded, r1_out, r2_out, compressed):
    '''
    (str, bool, str, str | None, bool) -> (_io.TextIOWrapper | None, _io.TextIOWrapper | None)
    
    Returns a tuple with fastq files opened for writing reads without matching
    patterns if keep_discarded is True or a tuple with None if False.
    Return a tuple with opened file and None if data is 'single'
    Fastq files are written in the directory of r1_out and r2_out and are named
    by appending '.non_matching_reads.RN.fastq.gz' to r1_out and r2_out
    
    Parameters
    ----------
    
    - data (str): Indicates if single or paired end sequencing data
    - keep_discarded (bool): Write reads without matching pattern to file if True
    - r1_out (str): Path to the output fastq 1 with reads re-headered with UMI sequence 
    - r2_out (str | None): Path to the output fastq 2 with reads re-headered with UMI sequence    
                           None for single end read sequences
    - compressed (bool): output_files r1_out and r2_out are compressed with gzip if True
    '''

    # initialize variables
    r1_discarded, r2_discarded = None, None

    # add suffix to file name 
    if compressed:
        suffix = '.discarded.{0}.fastq.gz'
    else:
        suffix = '.discarded.{0}.fastq'

    # open optional files for writing. same directory as output fastqs
    if keep_discarded:
        if data == 'paired':
            r1_discarded = _open_fastq_writing(_remove_fastq_extension(r1_out) + suffix.format('R1'), compressed)
            r2_discarded = _open_fastq_writing(_remove_fastq_extension(r2_out) + suffix.format('R2'), compressed)
        elif data == 'single':
            r1_discarded = _open_fastq_writing(_remove_fastq_extension(r1_out) + suffix.format('R1'), compressed)
    
    return r1_discarded, r2_discarded



def _remove_fastq_extension(fastq):
    '''
    (str) -> str
    
    Returns the name of the fastq file without extension or returns the file name
    if the extension is not one of: ".fastq.gz", ".fastq", ".fq.gz", ".fq", ".gz"
    
    Parameters
    ----------
    
    - fastq (str): Path to the fastq file
        
    Examples
    --------
    >>> _remove_fastq_extension('/myfolder/folder/folder/myfile.fastq.gz')
    '/myfolder/folder/folder/myfile'
    >>> _remove_fastq_extension('/myfolder/folder/folder/myfile.fastq.gz.fq.gz')
    '/myfolder/folder/folder/myfile'
    >>> _remove_fastq_extension('/myfolder/folder/folder/myfile.fastq.gz.fq.gz.fq')
    '/myfolder/folder/folder/myfile'
    >>> _remove_fastq_extension('myfile.R1.fastq.gz')
    'myfile.R1'
    >>> _remove_fastq_extension('myfile.R1.fastq.gz.trimmed.fq.gz')
    'myfile.R1.fastq.gz.trimmed'
    '''
    
    # copy the file name
    name = fastq
    # make a list of possible Fastq extensions
    extensions = ['.fastq.gz', '.fastq', '.fq.gz', '.fq', '.gz']
    while name[name.rfind('.'):] in extensions:
        name = name[: name.rfind('.')]
    
    return name



def _add_gzip_extension(r1_out, r2_out, compressed):
    '''
    (str, str | None, bool) -> (str, str)
    
    Returns file names r1_out and r2_out, if defined, with '.gz' extension if
    compressed is True and '.gz' not already in file name
    
    Parameters
    ----------
    - r1_out (str): Path to the output FASTQ 1 with reads re-headered with UMI sequence 
    - r2_out (str or None): Path to the output FASTQ 2 with reads re-headered with UMI sequence    
                            None for single end read sequences
    - compressed (bool): output fastqs are compressed with gzip if True
        
    Examples
    --------
    >>> _add_gzip_extension('myfile.R1.fastq.gz', 'myfile.R2.fastq.gz', False)
    ('myfile.R1.fastq.gz', 'myfile.R2.fastq.gz')
    >>> _add_gzip_extension('myfile.R1.fastq.gz', 'myfile.R2.fastq.gz', True)
    ('myfile.R1.fastq.gz', 'myfile.R2.fastq.gz')
    >>> _add_gzip_extension('myfile.R1.fastq.gz', None, True)
    ('myfile.R1.fastq.gz', None)
    >>> _add_gzip_extension('myfile.R1.fastq', None, True)
    ('myfile.R1.fastq.gz', None)
    >>> _add_gzip_extension('myfile.R1.fastq', 'myfile.R2.fastq', True)
    ('myfile.R1.fastq.gz', 'myfile.R2.fastq.gz')
    '''
    
    if compressed:
        if r1_out[-3:] != '.gz':
            r1_out += '.gz'
        if r2_out:
            if r2_out[-3:] != '.gz':
                r2_out += '.gz'
    return r1_out, r2_out


def _get_valid_barcodes(umilist):
    '''
    (str) -> list
    
    Returns a list of accepted barcodes present in column 1 of the umilist
    file. Any other columns, if present, are ignored.
    
    Parameters
    ----------
    - umilist (str): Path to file with accepted barcodes. Barcodes are expected
                       in the 1st column. Any other columns are ignored
    '''
    
    # open file with accepted barcodes
    infile = open(umilist)
    # create a list of barcodes
    barcodes = []
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split(' ')
            # barcodes expected in column 1. other columns ignored
            barcodes.append(line[0].strip())
    infile.close()
    return barcodes


def _write_discarded_reads(keep_discarded, discarded_fastqs, read):
    '''
    (bool, list, list) -> None
    
    Write reads without matching patterns to corresponding fastqs
    
    Parameters
    ----------
    - keep_discarded (bool): Write reads without matching patterns to file if True
    - discarded_fastqs (list): List of files opened for writing
    - read (list): List of reads
    '''
    
    if keep_discarded:
        for i in range(len(discarded_fastqs)):
            discarded_fastqs[i].write('\n'.join(list(map(lambda x: x.strip(), read[i]))) + '\n')


def _write_metrics(D, outputfile):
    '''
    (dict, str) -> None
    
    Writes data in D as a json
    
    Parameters
    ----------
    
    - D (dict): Data in the form of a dictionary of key, value pairs
    - outputfile (str): Path the output json file
    '''
    
    with open(outputfile, 'w') as newfile:
        json.dump(D, newfile, indent=4)


def extract_barcodes(r1_in, r1_out, pattern1, pattern2=None, inline_umi=True,
                     data='single', keep_extracted=True, keep_discarded=True,
                     r2_in=None, r2_out=None, r3_in=None, full_match=False,
                     separator='_', compressed=True, umilist=None):
    '''
    (list, str, str | None, str | None, bool, str, bool, bool, list | None, str | None, list | None, bool, str, bool, str | None) -> None

    Parameters
    ----------
    
    - r1_in (list): Path(s) to the input FASTQ 1 (compressed or not)
    - r1_out (str): Path to the output FASTQ 1 with reads re-headered with UMI sequence 
    - pattern1 (str or None): String sequence or regular expression used for matching and extracting UMis from reads in FASTQ 1.
                             The string sequence must look like NNNATCG or NNN. UMI nucleotides are labeled with "N".
                             Spacer nucleotides following Ns are used for matching UMIs but are discarded from reads    
                             None if UMIs are extracted only from FASTQ 2 for paired end sequences
    - pattern2 (str or None): String sequence or regular expression used for matching and extracting UMis from reads in FASTQ 2.
                             The string sequence must look like NNNATCG or NNN. UMI nucleotides are labeled with "N".
                             Spacer nucleotides following Ns are used for matching UMIs but are discarded from reads    
                             None if UMIs are extracted only from FASTQ 1 for paired end sequences
    - inline_umi (bool): True if UMIs are inline with reads and False otherwise
    - data (str): Indicates if single or paired end sequencing data
    - keep_extracted (bool): Write extracted sequences (UMIs and discarded sequences) to file if True
    - keep_discarded (bool): Write reads without matching pattern to file if True
    - r2_in (list or None): Path(s) to the input FASTQ 2 (compressed or not)
    - r2_out (str or None): Path to the output FASTQ 2 with reads re-headered with UMI sequence    
                           None for single end read sequences
    - r3_in (list or None): Path(s) to input FASTQ 3 for paired end sequences with non-inline UMIs 
    - full_match (bool): True if the regular expression needs to match the entire read sequence
    - separator (str): String separating the UMI sequence and part of the read header
    - compressed (bool): output fastqs are compressed with gzip if True
    - umilist (str or None): Path to file with accepted barcodes. Barcodes are expected
                               in the 1st column. Any other columns are ignored.
    '''
    
    # time function call
    start = time.time()

    # check that the number of input files is the same for paired end data
    _check_input_files(r1_in, r2_in, r3_in)
    
    # check input and output parameters
    _check_input_output(r1_in, r1_out, data, inline_umi, r2_in, r2_out, r3_in)
    # check pattern parameters 
    _check_pattern_options(pattern1, pattern2, data, inline_umi)

    # open outfiles for writing
    # add '.gz' if output fastqs are compressed and '.gz' not in file name
    r1_out, r2_out = _add_gzip_extension(r1_out, r2_out, compressed)    
    r1_writer = _open_fastq_writing(r1_out, compressed)
    r2_writer = _open_fastq_writing(r2_out, compressed) if data == 'paired' and r2_out else None
    
    # open optional files for writing. same directory as output fastqs
    # open files for writing reads without matching patterns
    r1_discarded, r2_discarded = _get_files_discarded_reads(data, keep_discarded, r1_out, r2_out, compressed)
    # open files for writing reads with extracted sequences (UMI and discarded sequences)
    r1_extracted, r2_extracted, r3_extracted = _get_files_extracted_reads(keep_extracted, data, inline_umi, pattern1, pattern2, r1_out, r2_out, compressed)
    
    # check that both patterns are either strings or regex
    _check_extraction_mode(pattern1, pattern2)
    # get pattern variables for each read 
    vals = list(zip(*list(map(lambda x: _get_read_patterns(x), [pattern1, pattern2]))))
    P = [pattern1, pattern2]
    seq_extract  = any(vals[0])
    UMIs = [vals[1][i] for i in range(len(P)) if P[i] is not None]
    spacers = [vals[2][i] for i in range(len(P)) if P[i] is not None]
    ps = [vals[3][i] for i in range(len(P)) if P[i] is not None]
    patterns = [i for i in P if i is not None]
    
    # make a list of files open for writing
    outfastqs = [i for i in [r1_writer, r2_writer] if i is not None]

    # make lists of optional files
    discarded_fastqs = [i for i in [r1_discarded, r2_discarded] if i is not None]
    extracted_fastqs = [i for i in [r1_extracted, r2_extracted, r3_extracted] if i is not None]
    
    # check if list of accepted barcodes provides
    if umilist:
        barcodes = _get_valid_barcodes(umilist)
    
    # count all reads and reads with matching and non-matching patterns
    Total, Matching, NonMatching = 0, 0, 0
    # track umi counts
    umi_counts = {}     
    
    # group input files
    file_groups = _group_input_files(r1_in, r2_in, r3_in)
    
    # loop over file groups, make list of opened files
    for group in file_groups:
        # open files for reading
        r1, r2, r3 = list(map(lambda x: _open_fastq(x) if x else None, group))
        # make a list of fastqs open for reading
        infastqs = [i for i in [r1, r2, r3] if i is not None]
         
        # create iterator with reads from each file
        Reads = zip(*map(lambda x: _get_read(x), infastqs))
    
        # loop over iterator with slices of 4 read lines from each file
        for read in Reads:
            # remove end of line from each read line
            read = _read_cleanup(read)
                
            # reset variable at each iteration. used to evaluate match
            umi = ''
            # check that input fastqs are in sync
            _check_fastq_sync([i[0] for i in read])
            # count total reads
            Total += 1
            # extract umis from reads
            if inline_umi:
                # extract UMI from read1 and/or read2 
                L = [_extract_umi_from_read(read[i], seq_extract, UMIs[i], spacers[i], ps[i], full_match) if patterns[i] else None for i in range(len(patterns))]     
                #L = [_extract_umi_from_read(read[i], seq_extract, UMIs[i], spacers[i], ps[i], full_match) for i in range(len(patterns)) if patterns[i]]     
            else:
                # UMIs are in fastq2 or fastq3 respectively for single and paired end data
                L = [_extract_umi_from_read(read[-1], seq_extract, UMIs[i], spacers[i], ps[i], full_match) if patterns[i] else None for i in range(len(patterns))]     
            # get umi sequences
            umi_sequences = [L[i][2] if L[i] else '' for i in range(len(L))]
            if all(map(lambda x: x is not None, L)) and all(map(lambda x: x != '', umi_sequences)):
                umi = ''.join(umi_sequences)
            
            # check if umi matched pattern
            if umi:
                # check if list of accepted barcodes
                # need to check if each umi+discarded seq is found in the list
                if umilist and all(map(lambda x: x in barcodes, umi_sequences)) == False:
                    # skip umis not listed
                    NonMatching += 1
                    # write non-matching reads to file if keep_discarded
                    _write_discarded_reads(keep_discarded, discarded_fastqs, read)
                else:
                    Matching +=1
                    # get read names, read, umi and extracted sequences and qualities for single and paired end
                    readnames = list(map(lambda x : _add_umi_to_readname(x, umi, separator), [read[i][0] for i in range(len(read))])) 
                    seqs, quals, umi_seqs, extracted_seqs, extracted_quals = zip(*L)
                
                    assert umi == ''.join(umi_seqs)
                    # update umi counter. keep track of UMI origin
                    umi_counts['.'.join(umi_seqs)] = umi_counts.get('.'.join(umi_seqs), 0) + 1
                                        
                    if inline_umi:
                                
                        if pattern1 and pattern2:
                            # paired end sequencing, umi extracted from each read
                            newreads = [[readnames[i], seqs[i], read[i][2], quals[i]] for i in range(len(read))]
                            # write extracted sequences to file(s)
                            if keep_extracted:
                                for i in range(len(extracted_fastqs)):
                                    extracted_fastqs[i].write('\n'.join([read[i][0], extracted_seqs[i], read[i][2], extracted_quals[i]]) + '\n')
                        elif pattern1:
                            # single or paired end sequencing, umi extracted from read1
                            newreads = [[readnames[0], seqs[0], read[0][2], quals[0]]]
                            if data == 'paired':
                                # no extraction from read2, append umi to read name and write read from input fastq2
                                newreads.append([readnames[1], read[1][1], read[1][2], read[1][3]])
                            if keep_extracted:
                                r1_extracted.write('\n'.join([read[0][0], extracted_seqs[0], read[0][2], extracted_quals[0]]) + '\n')
                        elif pattern2:
                            # paired end sequencing, umi extracted from read2
                            newreads = [list(map(lambda x: x.strip(), [readnames[0], read[0][1], read[0][2], read[0][3]]))]
                            newreads.append(list(map(lambda x: x.strip(), [readnames[1], seqs[0], read[1][2], quals[0]])))
                            if keep_extracted and r2_extracted:
                                r2_extracted.write('\n'.join([read[1][0], extracted_seqs[0], read[1][2], extracted_quals[0]]) + '\n')

                    else:
                        # single end: umi extracted from read 2. paired end: umi extracted from read 3
                        # keep read sequence and qualities
                        newreads = [[readnames[i], read[i][1], read[i][2], read[i][3]] for i in range(len(read) -1)]
                        if keep_extracted:
                            extracted_fastqs[-1].write('\n'.join([read[-1][0], extracted_seqs[0], read[-1][2], extracted_quals[0]]) + '\n')
                    
                    # write new reads to output fastq
                    for i in range(len(outfastqs)):
                        outfastqs[i].write('\n'.join(newreads[i]) +'\n')
            else:
                NonMatching += 1
                # write non-matching reads to file if keep_discarded
                _write_discarded_reads(keep_discarded, discarded_fastqs, read)
    
        # close all input fastqs
        for i in infastqs:
            i.close()
        
    # close all output fastqs
    for i in outfastqs + discarded_fastqs + extracted_fastqs:
        i.close()
    
    # save metrics to files
    d = {'total reads/pairs': Total, 'reads/pairs with matching pattern': Matching,
         'discarded reads/pairs': NonMatching, 'pattern1': pattern1, 'pattern2': pattern2,
         'umi-list file': umilist}
    # get prefix from r1_out
    prefix = _remove_fastq_extension(os.path.basename(r1_out))
    _write_metrics(d, os.path.join(os.path.dirname(r1_out), '{0}_extraction_metrics.json'.format(prefix)))
    _write_metrics(umi_counts, os.path.join(os.path.dirname(r1_out), '{0}_UMI_counts.json'.format(prefix)))
    
    print('total reads/pairs:', Total)
    print('reads/pairs with matching pattern:', Matching)
    print('discarded reads/pairs:', NonMatching)
    if pattern1:
        print('pattern1:', pattern1)
    if pattern2:
        print('pattern2', pattern2)
    if umilist:
        print('umi-list file:', umilist)
    
    # record time after call
    end = time.time()
    print('Extracted UMIs in {0} seconds'.format(round(end - start, 3)))



    
def main():
        
    # create parser
    parser = argparse.ArgumentParser(prog='barcodex.py', description='A package for extracting Unique Molecular Identifiers (UMIs) from single or paired read sequencing data')
    subparsers = parser.add_subparsers()
       		
    # extract commands
    e_parser = subparsers.add_parser('extract', help="Extract UMIs from read sequences")
    e_parser.add_argument('--r1_in', dest='r1_in', nargs='*', help='Path to input FASTQ 1', required=True)
    e_parser.add_argument('--r1_out', dest='r1_out', help='Path to output FASTQ 1', required=True)
    e_parser.add_argument('--pattern1', dest='pattern1', help='Barcode string of regex for extracting UMIs in read 1')
    e_parser.add_argument('--pattern2', dest='pattern2', help='Barcode string of regex for extracting UMIs in read 2')
    e_parser.add_argument('--inline', dest='inline_umi', action='store_true', help='UMIs inline with reads or not. True if activated')
    e_parser.add_argument('--data', dest='data', choices=['single', 'paired'], default='single', help='Paired or single end sequencing')
    e_parser.add_argument('--r2_in', dest='r2_in', nargs='*', help='Path to input FASTQ 2. Fastq 2 for paired end sequencing with inline UMIs. Fastq with UMIs for single end sequencing with UMIs not in line')
    e_parser.add_argument('--r2_out', dest='r2_out', help='Path to output FASTQ 2')
    e_parser.add_argument('--r3_in', dest='r3_in', nargs='*', help='Path to input FASTQ 3. Fastq with UMIs for paired end sequencing with UMIs not in line')
    e_parser.add_argument('--separator', dest='separator', default='_', help='String separating the UMI sequence in the read name')
    e_parser.add_argument('--keep_extracted', dest='keep_extracted', action='store_true', help='Output the extracted UMIs and potentially discarded sequences from reads in separate fastqs. True if activated')
    e_parser.add_argument('--keep_discarded', dest='keep_discarded', action='store_true', help='Output reads with non-matching patterns to separate fastqs. True if activated')
    e_parser.add_argument('--full_match', dest='full_match', action='store_true', help='Requires the regex pattern to match the entire read sequence. True if activated')
    e_parser.add_argument('--compressed', dest='compressed', action='store_true', help='Compress output fastqs with gzip. True if activated')
    e_parser.add_argument('--umilist', dest='umilist', help='Path to file with valid UMIs (1st column)')
        
    args = parser.parse_args()
    
    try:
        extract_barcodes(args.r1_in, args.r1_out, pattern1=args.pattern1, pattern2=args.pattern2,
                     inline_umi=args.inline_umi, data=args.data, keep_extracted=args.keep_extracted, keep_discarded=args.keep_discarded,
                     r2_in=args.r2_in, r2_out=args.r2_out, r3_in=args.r3_in, full_match=args.full_match, compressed=args.compressed, umilist=args.umilist)
    except AttributeError as e:
        print('#############\n')
        print('AttributeError: {0}\n'.format(e))
        print('#############\n\n')
        print(parser.format_help())
 
