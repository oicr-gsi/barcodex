# -*- coding: utf-8 -*-
"""
Created on Tue May 12 14:23:08 2020

@author: rjovelin
"""

import gzip
from itertools import zip_longest
import regex
import argparse
import json
import time



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
        
    
def _check_input_files(L1, L2, L3=None):
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


def _group_input_files(L1, L2, L3=None):
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
    

def _check_pattern_options(r1_in, pattern1, r2_in=None, pattern2=None):
    '''
    (list, str | None, list | None, str | None) -> None
    
    Raise ValueError if pattern options are incompatible with single or paired
    sequencing data
    
    Parameters
    ----------
    - r1_in (list): Paths to FASTQs 1
    - pattern1 (str or None): String sequence or regular expression used for matching and extracting UMis from reads in FASTQ 1.
                              None if UMIs are extracted only from FASTQ 2 
    - r2_in (list or None): Paths to FASTQs 2 for paired end. None for single end 
    - pattern2 (str or None): String sequence or regular expression used for matching and extracting UMis from reads in FASTQ 2.
                              None if UMIs are extracted only from FASTQ 1.
        
    Examples
    --------    
    >>> _check_pattern_options(['file1.fatsq'], 'NNNATCG', None, None)
    >>> _check_pattern_options(['file1.fatsq'], 'NNNATCG', ['file2.fastq'], None)
    >>> _check_pattern_options(['file1.fatsq'], None, ['file2.fastq'], 'NNNATCG')
    >>> _check_pattern_options(['file1.fatsq'], None, ['file2.fastq'], None)
    ValueError: Expecting paired end sequences. At least 1 pattern is required
    >>> _check_pattern_options(['file1.fatsq'], None, None, None)
    ValueError: Expecting single end sequences. Pattern1 required, pattern2 not needed
    '''

    if r1_in and r2_in is None:
        # single end. expecting pattern but not pattern2
        if pattern2 or pattern1 is None:
            raise ValueError('Expecting single end sequences. Pattern1 required, pattern2 not needed')
    elif r1_in and r2_in:
        # paired end. at least 1 pattern must be used
        if pattern1 is None and pattern2 is None:
            raise ValueError('Expecting paired end sequences. At least 1 pattern is required')
    
         
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
    - a boolean indicating whether UMI is extracted with a string pattern (True) or regex (False)
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


def _open_fastq_writing(output_file):
    '''
    (str) -> _io.TextIOWrapper
    
    Returns a file handler for writing to output_file. The output_file is compressed
    with gzip (highest compression, level 9)
    
    Parameters
    ----------
    
    - output_file (str): Path to the output_file
    '''
    
    newfile = gzip.GzipFile(filename=None, mode="w", fileobj=open(output_file, 'wb'), mtime=0)
    return newfile


def _get_valid_barcodes(umilist):
    '''
    (str) -> list
    
    Returns a list of accepted barcodes.
    
    Parameters
    ----------
    - umilist (str): Path to file with accepted barcodes. 
   '''
    
    # set lists for barcoes found in read1 and read 2
    barcode1 = []
    barcode2 = []
    valid_umis = []    
    
    # open file with accepted barcodes
    infile = open(umilist)
    for line in infile:
        line = line.rstrip()
        if line != '':
            parts = line.split()
            if len(parts) == 1:
                # umi expected in read1 and read2
                barcode1.append(parts[0])
                barcode2.append(parts[0])
            elif len(parts) > 1:
                # read is indicated by 1 or 2 or both
                for p in parts[1:]:
                    if p == "1":
                        # barcode expected in read1
                        barcode1.append(parts[0])
                    elif p == "2":
                        # barcode expected in read2
                        barcode2.append(parts[0])

    if len(barcode1) == 0:
        valid_umis = barcode2
    elif len(barcode2) == 0:
        valid_umis = barcode1
    else:
        for i in barcode1:
            for j in barcode2:
                valid_umis.append(i + '.' + j)
    return valid_umis


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




def _get_extraction_metrics(Total, Matching, NonMatching, Rejected, pattern1=None, pattern2=None, umilist=None):
    '''
    (int, int, int, int, str | None, str | None, str | None) -> dict
    
    Parameters
    ----------
    - Total (int): Total reads/pairs
    - Matching (int): Reads/pairs with matching pattern
    - NonMatching (int): Discarded reads/pairs
    - Rejected (int): Discarded reads/pairs due to unknown UMI
    - pattern1 (str or None): Pattern used to extract reads in FASTQ1
    - pattern2 (str or None): Pattern used to extract reads in FASTQ2
    - umilist (str or None): Path to file with valid UMIs
    
    Returns a dictionary with extraction metrics and print metrics
    '''
    
    # save metrics to files
    d = {'total reads/pairs': Total, 'reads/pairs with matching pattern': Matching,
         'discarded reads/pairs': NonMatching, 'discarded reads/pairs due to unknown UMI': Rejected}
    
    print('total reads/pairs:', Total)
    print('reads/pairs with matching pattern:', Matching)
    print('discarded reads/pairs:', NonMatching)
    print('discarded reads/pairs due to unknown UMI:', Rejected)
    if pattern1:
        print('pattern1:', pattern1)
        d['pattern1'] = pattern1
    if pattern2:
        print('pattern2', pattern2)
        d['pattern2'] = pattern2
    if umilist:
        print('umi-list file:', umilist)
        d['umi-list file'] = umilist
       
    return d

    
def extract_barcodes_inline(r1_in, pattern1, prefix, pattern2=None, r2_in=None,
                            full_match=False, separator='_', umilist=None):
    '''
    (list, str | None, str, str | None, str | None, bool, str, str | None) -> None

    Parameters
    ----------
    - r1_in (list): Path(s) to the input FASTQ 1
    - pattern1 (str or None): String sequence or regular expression used for matching and extracting UMis from reads in FASTQ 1.
                             The string sequence must look like NNNATCG or NNN. UMI nucleotides are labeled with "N".
                             Spacer nucleotides following Ns are used for matching UMIs but are discarded from reads    
                             None if UMIs are extracted only from FASTQ 2 for paired end sequences
    - prefix (str): Specifies the start of the output files and stats json files
    - pattern2 (str or None): String sequence or regular expression used for matching and extracting UMis from reads in FASTQ 2.
                             The string sequence must look like NNNATCG or NNN. UMI nucleotides are labeled with "N".
                             Spacer nucleotides following Ns are used for matching UMIs but are discarded from reads    
                             None if UMIs are extracted only from FASTQ 1 for paired end sequences
    - r2_in (list or None): Path(s) to the input FASTQ 2
    - full_match (bool): True if the regular expression needs to match the entire read sequence
    - separator (str): String separating the UMI sequence and part of the read header
    - umilist (str or None): Path to file with accepted barcodes
    '''
    
    # time function call
    start = time.time()

    # check that the number of input files is the same for paired end data
    _check_input_files(r1_in, r2_in)
    
    # check pattern parameters 
    _check_pattern_options(r1_in, pattern1, r2_in=r2_in, pattern2=pattern2)

    # open outfiles for writing
    r1_writer = _open_fastq_writing(prefix + '_R1.fastq.gz')
    r2_writer = _open_fastq_writing(prefix + '_R2.fastq.gz') if r2_in else None
    
    # open optional files for writing. same directory as output fastqs
    # open files for writing reads without matching patterns
    r1_discarded, r2_discarded = None, None
    r1_discarded = _open_fastq_writing(prefix + '_discarded.R1.fastq.gz')
    if r2_in:
        # paired data
        r2_discarded = _open_fastq_writing(prefix + '_discarded.R2.fastq.gz')
    
    # open files for writing reads with extracted sequences (UMI and discarded sequences)
    r1_extracted, r2_extracted = None, None
    if pattern1 is not None:
        r1_extracted = _open_fastq_writing(prefix +  '_extracted.R1.fastq.gz')
    if pattern2 is not None:
            r2_extracted = _open_fastq_writing(prefix +  '_extracted.R2.fastq.gz')
    
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
    extracted_fastqs = [i for i in [r1_extracted, r2_extracted] if i is not None]
    
    # check if list of accepted barcodes provides
    if umilist:
        barcodes = _get_valid_barcodes(umilist)
    
    # count all reads and reads with matching and non-matching patterns
    Total, Matching, NonMatching, Rejected = 0, 0, 0, 0
    # track umi counts
    umi_counts = {}     
    
    # group input files
    file_groups = _group_input_files(r1_in, r2_in)
    
    # loop over file groups, make list of opened files
    for group in file_groups:
        # open files for reading
        r1, r2, _ = list(map(lambda x: gzip.open(x, 'rt') if x else None, group))
        # make a list of fastqs open for reading
        infastqs = [i for i in [r1, r2] if i is not None]
         
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
            # extract UMI from read1 and/or read2 
            L = [_extract_umi_from_read(read[i], seq_extract, UMIs[i], spacers[i], ps[i], full_match) if patterns[i] else None for i in range(len(patterns))]     
            #L = [_extract_umi_from_read(read[i], seq_extract, UMIs[i], spacers[i], ps[i], full_match) for i in range(len(patterns)) if patterns[i]]     
                
            # get umi sequences
            umi_sequences = [L[i][2] if L[i] else '' for i in range(len(L))]
            if all(map(lambda x: x is not None, L)) and all(map(lambda x: x != '', umi_sequences)):
                # use dot to join UMIs from each read
                umi = '.'.join(umi_sequences)
            
            # check if umi matched pattern
            if umi:
                # check if list of accepted barcodes
                # check if concatenated umi in list of valid barcodes
                if umilist and umi not in barcodes:
                    # skip umis not listed
                    Rejected += 1
                    # write non-matching reads to file
                    for i in range(len(discarded_fastqs)):
                        discarded_fastqs[i].write(str.encode('\n'.join(list(map(lambda x: x.strip(), read[i]))) + '\n'))
                else:
                    Matching +=1
                    # get read names, read, umi and extracted sequences and qualities for single and paired end
                    readnames = list(map(lambda x : _add_umi_to_readname(x, umi, separator), [read[i][0] for i in range(len(read))])) 
                    seqs, quals, umi_seqs, extracted_seqs, extracted_quals = zip(*L)
                
                    assert umi == '.'.join(umi_seqs)
                    # update umi counter. keep track of UMI origin
                    umi_counts['.'.join(umi_seqs)] = umi_counts.get('.'.join(umi_seqs), 0) + 1
                                        
                    if pattern1 and pattern2:
                        # paired end sequencing, umi extracted from each read
                        newreads = [[readnames[i], seqs[i], read[i][2], quals[i]] for i in range(len(read))]
                        # write extracted sequences to file(s)
                        for i in range(len(extracted_fastqs)):
                            extracted_fastqs[i].write(str.encode('\n'.join([read[i][0], extracted_seqs[i], read[i][2], extracted_quals[i]]) + '\n'))
                    elif pattern1:
                        # single or paired end sequencing, umi extracted from read1
                        newreads = [[readnames[0], seqs[0], read[0][2], quals[0]]]
                        if r2_in:
                            # paired data. no extraction from read2, append umi to read name and write read from input fastq2
                            newreads.append([readnames[1], read[1][1], read[1][2], read[1][3]])
                        r1_extracted.write(str.encode('\n'.join([read[0][0], extracted_seqs[0], read[0][2], extracted_quals[0]]) + '\n'))
                    elif pattern2:
                        # paired end sequencing, umi extracted from read2
                        newreads = [list(map(lambda x: x.strip(), [readnames[0], read[0][1], read[0][2], read[0][3]]))]
                        newreads.append(list(map(lambda x: x.strip(), [readnames[1], seqs[0], read[1][2], quals[0]])))
                        if r2_extracted:
                            r2_extracted.write(str.encode('\n'.join([read[1][0], extracted_seqs[0], read[1][2], extracted_quals[0]]) + '\n'))
                    # write new reads to output fastq
                    for i in range(len(outfastqs)):
                        outfastqs[i].write(str.encode('\n'.join(newreads[i]) +'\n'))
            else:
                NonMatching += 1
                # write non-matching reads to file
                for i in range(len(discarded_fastqs)):
                    discarded_fastqs[i].write(str.encode('\n'.join(list(map(lambda x: x.strip(), read[i]))) + '\n'))

        # close all input fastqs
        for i in infastqs:
            i.close()
        
    # close all output fastqs
    for i in outfastqs + discarded_fastqs + extracted_fastqs:
        i.close()

    # get extraction metrics
    d = _get_extraction_metrics(Total, Matching, NonMatching, Rejected, pattern1=pattern1, pattern2=pattern2, umilist=umilist)
    # save metrics to files
    _write_metrics(d, prefix + '_extraction_metrics.json')
    _write_metrics(umi_counts, prefix + '_UMI_counts.json')
    
    # record time after call
    end = time.time()
    print('Extracted UMIs in {0} seconds'.format(round(end - start, 3)))


def extract_barcodes_separate(r1_in, pattern1, prefix, ru_in, r2_in=None, full_match=False, separator='_', umilist=None):
    '''
    (list, str, str, list, list | None, bool, str, str | None) -> None

    Parameters
    ----------

    - r1_in (list): Path(s) to the input FASTQ 1
    - pattern1 (str): String sequence or regular expression used for matching and extracting UMis
                      from reads in FASTQ 2 (single end )or FASTQ2 (paired end). 
                      The string sequence must look like NNNATCG or NNN. UMI nucleotides are labeled with "N".
                      Spacer nucleotides following Ns are used for matching UMIs but are discarded from reads
    - prefix (str): Specifies the start of the output files and stats json files
    - ru_in (list): Path(s) to input FASTQ containing UMis
    - r2_in (list or None): Path(s) to input FASTQ 2 for paired  
    - full_match (bool): True if the regular expression needs to match the entire read sequence
    - separator (str): String separating the UMI sequence and part of the read header
    - umilist (str or None): Path to file with accepted barcodes
    '''
    
    # time function call
    start = time.time()

    # check that the number of input files is the same for paired end data
    _check_input_files(r1_in, ru_in, r2_in)
    
    # open outfiles for writing
    r1_writer = _open_fastq_writing(prefix + '_R1.fastq.gz')
    r2_writer = _open_fastq_writing(prefix + '_R2.fastq.gz') if r2_in else None
    
    # open optional files for writing. same directory as output fastqs
    # open files for writing reads without matching patterns
    r1_discarded, r2_discarded = None, None
    r1_discarded = _open_fastq_writing(prefix + '_discarded.R1.fastq.gz')
    if r2_in:
        # data is paired. open file to discard read 2
        r2_discarded = _open_fastq_writing(prefix + '_discarded.R2.fastq.gz')
    # open files for writing reads with extracted sequences (UMI and discarded sequences)
    ru_extracted = _open_fastq_writing(prefix + '_extracted.fastq.gz')
        
    # get pattern variables for each read 
    seq_extract, UMI, spacer, ps = _get_read_patterns(pattern1)
        
    # make a list of files open for writing
    outfastqs = [i for i in [r1_writer, r2_writer] if i is not None]

    # make lists of optional files
    discarded_fastqs = [i for i in [r1_discarded, r2_discarded] if i is not None]
        
    # check if list of accepted barcodes provides
    if umilist:
        barcodes = _get_valid_barcodes(umilist)
    
    # count all reads and reads with matching and non-matching patterns
    Total, Matching, NonMatching, Rejected = 0, 0, 0, 0
    # track umi counts
    umi_counts = {}     
    
    # group input files
    file_groups = _group_input_files(r1_in, ru_in, r2_in)
    
    # loop over file groups, make list of opened files
    for group in file_groups:
        # open files for reading
        r1, ru, r2 = list(map(lambda x: gzip.open(x, 'rt') if x else None, group))
        # make a list of fastqs open for reading
        infastqs = [i for i in [r1, ru, r2] if i is not None]
         
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
            # UMIs are in fastq2 or fastq3 respectively for single and paired end data
            L = _extract_umi_from_read(read[-1], seq_extract, UMI, spacer, ps, full_match)
            # get umi sequences
            umi = L[2]
            
            # check if umi matched pattern
            if umi:
                # check if list of accepted barcodes
                # check if concatenated umi in list of valid barcodes
                if umilist and umi not in barcodes:
                    # skip umis not listed
                    Rejected += 1
                    # write non-matching reads to file
                    for i in range(len(discarded_fastqs)):
                        discarded_fastqs[i].write(str.encode('\n'.join(list(map(lambda x: x.strip(), read[i]))) + '\n'))
                else:
                    Matching +=1
                    # get read names, read, umi and extracted sequences and qualities for single and paired end
                    readnames = list(map(lambda x : _add_umi_to_readname(x, umi, separator), [read[i][0] for i in range(len(read))])) 
                    seqs, quals, umi_seqs, extracted_seqs, extracted_quals = L
                    assert umi == umi_seqs
                    # update umi counter. keep track of UMI origin
                    umi_counts[umi] = umi_counts.get(umi, 0) + 1
                                        
                    # single end: umi extracted from read 2. paired end: umi extracted from read 3
                    # keep read sequence and qualities
                    newreads = [[readnames[i], read[i][1], read[i][2], read[i][3]] for i in range(len(read) -1)]
                    ru_extracted.write(str.encode('\n'.join([read[-1][0], extracted_seqs, read[-1][2], extracted_quals]) + '\n'))
                    
                    # write new reads to output fastq
                    for i in range(len(outfastqs)):
                        outfastqs[i].write(str.encode('\n'.join(newreads[i]) +'\n'))
            else:
                NonMatching += 1
                # write non-matching reads to file
                for i in range(len(discarded_fastqs)):
                    discarded_fastqs[i].write(str.encode('\n'.join(list(map(lambda x: x.strip(), read[i]))) + '\n'))

    
        # close all input fastqs
        for i in infastqs:
            i.close()
        
    # close all output fastqs
    for i in outfastqs + discarded_fastqs + [ru_extracted]:
        i.close()
    
    # get extraction metrics
    d = _get_extraction_metrics(Total, Matching, NonMatching, Rejected, pattern1=pattern1, pattern2=None, umilist=umilist)
    # save metrics to files
    _write_metrics(d, prefix + '_extraction_metrics.json')
    _write_metrics(umi_counts, prefix + '_UMI_counts.json')
    
    # record time after call
    end = time.time()
    print('Extracted UMIs in {0} seconds'.format(round(end - start, 3)))

    
def main():
    # create parser
    parser = argparse.ArgumentParser(prog='barcodex.py', description='A package for extracting Unique Molecular Identifiers (UMIs) from single or paired read sequencing data')
    parser.add_argument('--umilist', dest='umilist', help='Path to file with valid UMIs (1st column)')
    parser.add_argument('--prefix', dest='prefix', help='Specifies the start of the output files and stats json files', required=True)
    parser.add_argument('--separator', dest='separator', default='_', help='String separating the UMI sequence in the read name')
    subparsers = parser.add_subparsers()
       		
    # inline UMIs 
    i_parser = subparsers.add_parser('inline', help="Extract UMIs located in read sequences")
    i_parser.add_argument('--r1_in', dest='r1_in', nargs='*', help='Path to input FASTQ 1', required=True)
    i_parser.add_argument('--pattern1', dest='pattern1', help='Barcode string of regex for extracting UMIs in read 1')
    i_parser.add_argument('--pattern2', dest='pattern2', help='Barcode string of regex for extracting UMIs in read 2')
    i_parser.add_argument('--r2_in', dest='r2_in', nargs='*', help='Path to input FASTQ 2. Fastq 2 for paired end sequencing with inline UMIs. Fastq with UMIs for single end sequencing with UMIs not in line')
    i_parser.add_argument('--full_match', dest='full_match', action='store_true', help='Requires the regex pattern to match the entire read sequence. True if activated')
       
    # UMI in separate file
    s_parser = subparsers.add_parser('separate', help="Extract UMIs located in separate file")
    s_parser.add_argument('--r1_in', dest='r1_in', nargs='*', help='Path to input FASTQ 1', required=True)
    s_parser.add_argument('--pattern1', dest='pattern1', help='Barcode string of regex for extracting UMIs in read 1')
    s_parser.add_argument('--ru_in', dest='ru_in', nargs='*', help='Path to input FASTQ containing UMIs')
    s_parser.add_argument('--r2_in', dest='r2_in', nargs='*', help='Path to input FASTQ 2. Fastq with UMIs for paired end sequencing with UMIs not in line')
    s_parser.add_argument('--full_match', dest='full_match', action='store_true', help='Requires the regex pattern to match the entire read sequence. True if activated')
    
    args = parser.parse_args()
    
    if args.subparser_name == 'inline':
        try:
            extract_barcodes_inline(args.r1_in, pattern1=args.pattern1, prefix=args.prefix, pattern2=args.pattern2, 
                                    r2_in=args.r2_in, separator=args.separator, umilist=args.umilist)
        except AttributeError as e:
            print('#############\n')
            print('AttributeError: {0}\n'.format(e))
            print('#############\n\n')
            print(parser.format_help())
    elif args.subparser_name == 'separate':
        try:
            extract_barcodes_separate(args.r1_in, pattern1=args.pattern1, prefix=args.prefix, ru_in=args.ru_in, r2_in=args.r2_in, full_match=args.full_match, separator=args.separator, umilist=args.umilist)
        except AttributeError as e:
            print('#############\n')
            print('AttributeError: {0}\n'.format(e))
            print('#############\n\n')
            print(parser.format_help())
    elif args.subparser_name is None:
        print(parser.format_help())
    