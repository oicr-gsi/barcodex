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

    :param pattern: A string sequence with UMI nucleotides as Ns and optional
                    spacer sequence. Pattern must look like NNNNATCG or NNNN

    Returns a list with UMI sequence and eventually the spacer sequence
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


def _check_extraction_mode(pattern, pattern2):
    '''
    (str | None, str | None) -> None
    
    Raise a ValueError if pattern and pattern2 are not both a string sequence or a regex
    
    Parameters
    ----------
    - pattern (str or None): String sequence or regular expression used for matching and extracting UMis from reads in FASTQ 1.
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
    
    if pattern and pattern2:
        if _is_pattern_sequence(pattern) == True and _is_pattern_sequence(pattern2) == False:
            raise ValueError('Both patterns must be either string sequences or regex')
        elif _is_pattern_sequence(pattern) == False and _is_pattern_sequence(pattern2) == True:
            raise ValueError('Both patterns must be either string sequences or regex')


def _get_umi_spacer(pattern):
    '''
    (str) -> (str, str)
    
    :param pattern: A string sequence with UMI nucleotides as Ns and optional
                    spacer sequence. Pattern must look like NNNNATCG or NNNN

    Return the UMI and spacer sequences in pattern 
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
    the umi sequence, the extracted read sequence and qualities
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
    (list, _regex.Pattern, bool)
    
    '''
    
    # initialize variables
    seq, qual, umi_seq, extracted_seq, extracted_qual = '', '', '', '', ''
        
    # look for a match in read sequence
    if full_match == False:
        # scan through the string looking for a match
        m = p.search(read[1])
    elif full_match == True:
        # match if the whole string matches pattern 
        m = full_match(read[1])
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
        removed = []
        for i in umi_pos + discard_pos:
            removed.extend(list(range(i[0], i[1])))
        remaining = sorted([i for i in range(len(read[1])) if i not in removed])
        
        # get extracted sequence and qualities
        extracted_seq, extracted_qual = '', ''
        for i in extracted_pos:
            extracted_seq += read[1][i[0]: i[1]]
            extracted_qual += read[3][i[0]: i[1]]
        # get read seq and qual after extraction
        seq, qual = '', ''
        for i in remaining:
            seq += read[1][i]
            qual += read[3][i]
        
    return seq, qual, umi_seq, extracted_seq, extracted_qual


def _get_read(fastq_file):
    """
    (_io.TextIOWrapper) -- > itertools.zip_longest
    :param fastq_file: a fastq file open for reading in plain text mode
    
    Returns an iterator slicing the fastq into 4-line reads.
    Each element of the iterator is a tuple containing read information
    """
    args = [iter(fastq_file)] * 4
    return zip_longest(*args, fillvalue=None)



def _check_fastq_sync(L):
    '''
    (list) -> None
    
    
    '''
    
    readnames = []
    for i in L:
        readnames.append(i.split(' ')[0])
    if len(set(readnames)) > 1:
        raise ValueError('Fastqs are not synced')
        
    
def _check_input_output(r1_in, r1_out, data='single', inline_umi=True,
                        r2_in=None, r2_out=None, r3_in=None):
    '''
    
    
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
            # requires pattern, pattern2 not needed
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


def _check_pattern_options(pattern, pattern2=None, data='single', inline_umi=True):
    '''
    
    
    '''
    
    
    if pattern is None and pattern2 is None:
        raise ValueError('At least 1 pattern is required')
    if not inline_umi:
        # pattern2 not needed
        if pattern2 is not None:
            raise ValueError('Expecting paired end sequences with out of read UMIs. Pattern2 not needed')
 
    
def _convert_arg_to_bool(argument):
    '''
    (str) -> bool
    
    "param argument: Argument of a parameter that should be a boolean
    
    Return the argument as a boolean 
    '''
    
    if isinstance(argument, bool):
       return argument
    elif argument.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif argument.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise ValueError('ERR: {0} is expected to be a boolean'.format(argument))
        
    
     
def _extract_umi_from_read(read, seq_extract, UMI, spacer, p, full_match):    
    '''
    (list, bool, str | None, str | None, _regex.Pattern | None, bool) -> (str, str, str, str, str)
    
    
    
    '''
    
    if seq_extract == True:
        # extraction using string sequence. assumes UMi starts at begining of read        
        L = _extract_from_sequence(read, UMI, spacer)
    else:
        L = _extract_from_regex(read, p, full_match)
    return L
    


def _get_read_patterns(pattern):
    '''
    
    
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





    

def extract_barcodes(r1_in, r1_out, pattern, pattern2=None, inline_umi=True,
                     data='single', keep_extracted=True, keep_discarded=True,
                     r2_in=None, r2_out=None, r3_in=None, full_match=False, separator='_'):
    """





    """
    
    # to do/things to consider
    # remove extension before appending discarded or extracted extension
    # append extracted sequences to read name instead of fastq
    # option to output to fastq or gzip compressed fastq.gz
    # option of compression level
    # run as module and script
    # handle umis not in line with read
    # use whitelist
    # add docstrings
    # add test cases
    
    
    
    # check input and output parameters
    _check_input_output(r1_in, r1_out, data, inline_umi, r2_in, r2_out, r3_in)
    # check pattern parameters 
    _check_pattern_options(pattern, pattern2, data, inline_umi)

    # open files for reading
    r1, r2, r3 = list(map(lambda x: _open_fastq(x) if x else None, [r1_in, r2_in, r3_in]))
    
    # open outfiles for writing
    r1_writer = gzip.open(r1_out, 'wt')
    if data == 'paired' and r2_out:
        r2_writer = gzip.open(r2_out, "wt")
    else:
        r2_writer = None
    
    # open optional files for writing. same directory as output fastqs
    if keep_discarded:
        # initialize variables
        r1_discarded, r2_discarded = None, None
        if data == 'paired':
            r1_discarded = gzip.open(r1_out + '.non_matching_reads.R1.fastq.gz', 'wt')
            r2_discarded = gzip.open(r2_out + '.non_matching_reads.R2.fastq.gz', 'wt')
        elif data == 'single':
            r1_discarded = gzip.open(r1_out + '.non_matching_reads.R1.fastq.gz', 'wt')
    else:
        r1_discarded, r2_discarded = None, None
    if keep_extracted:
        # initialize variables
        r1_extracted, r2_extracted, r3_extracted = None, None, None
        if inline_umi:
            if pattern is not None:
                r1_extracted = gzip.open(r1_out + '.extracted_sequences.R1.fastq.gz', 'wt')
            if pattern2 is not None:
                r2_extracted = gzip.open(r2_out + '.extracted_sequences.R2.fastq.gz', 'wt')
        else:
            if data == 'paired':
                outdir = os.path.direname(r1_out)
                filename = os.path.basename(r3_in)
                outfile = os.path.join(outdir, filename + '.umi_sequences.R3.fastq.gz')
                r3_extracted = gzip.open(outfile, 'wt')
            elif data == 'single':
                outdir = os.path.direname(r1_out)
                filename = os.path.basename(r2_in)
                outfile = os.path.join(outdir, filename + '.umi_sequences.R2.fastq.gz')
                r2_extracted = gzip.open(outfile, 'wt')
    else:
        r1_extracted, r2_extracted, r3_extracted = None, None, None
        
        
        
    # check that both patterns are either strings or regex
    _check_extraction_mode(pattern, pattern2)
    # get pattern variables for each read 
    vals = list(zip(*list(map(lambda x: _get_read_patterns(x), [pattern, pattern2]))))
    P = [pattern, pattern2]
    seq_extract  = any(vals[0])
    UMIs = [vals[1][i] for i in range(len(P)) if P[i] is not None]
    spacers = [vals[2][i] for i in range(len(P)) if P[i] is not None]
    ps = [vals[3][i] for i in range(len(P)) if P[i] is not None]
    patterns = [i for i in P if i is not None]
    
    # make a list of fastqs open for reading
    infastqs = [i for i in [r1, r2, r3] if i is not None]
    
    # make a list of files open for writing
    outfastqs = [i for i in [r1_writer, r2_writer] if i is not None]

    # make lists of optional files
    discarded_fastqs = [i for i in [r1_discarded, r2_discarded] if i is not None]
    extracted_fastqs = [i for i in [r1_extracted, r2_extracted, r3_extracted] if i is not None]
    
    # create iterator with reads from each file
    Reads = zip(*map(lambda x: _get_read(x), infastqs))
    
    # count all reads and reads with matching and non-matching patterns
    Total, Matching, NonMatching = 0, 0, 0
         
    # loop over iterator with slices of 4 read lines from each file
    for read in Reads:
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
            # get umi sequences
            umi_sequences = [L[i][2] if L[i] else '' for i in range(len(L))]
            if all(map(lambda x: x is not None, L)) and all(map(lambda x: x != '', umi_sequences)):
                umi = ''.join(umi_sequences)
            
        # check if umi matched pattern
        if umi:
            Matching +=1
            if inline_umi:
                # get read names, read, umi and extracted sequences and qualities for single and paired end
                readnames = list(map(lambda x : _add_umi_to_readname(x, umi, separator), [read[i][0] for i in range(len(read))])) 
                seqs, quals, umi_seqs, extracted_seqs, extracted_quals = zip(*L)
                
                assert umi == ''.join(umi_seqs)
                
                
                if pattern and pattern2:
                    # paired end sequencing, umi extracted from each read
                    newreads = [list(map(lambda x: x.strip(), [readnames[i], seqs[i], read[i][2], quals[i]])) for i in range(len(read))]
                    # write extracted sequences to file(s)
                    if keep_extracted:
                        for i in range(len(extracted_fastqs)):
                            extracted_fastqs[i].write('\n'.join(list(map(lambda x: x.strip(), [read[i][0], extracted_seqs[i], read[i][2], extracted_quals[i]]))))
                elif pattern:
                    # single or paired end sequencing, umi extracted from read1
                    newreads = [list(map(lambda x: x.strip(), [readnames[0], seqs[0], read[0][2], quals[0]]))]
                    if data == 'paired':
                        # no extraction from read2, append umi to read name and write read from input fastq2
                        newreads.append(list(map(lambda x: x.strip(), [readnames[1], read[1][1], read[1][2], read[1][3]])))
                    if keep_extracted:
                        r1_extracted.write('\n'.join(list(map(lambda x: x.strip(), [read[0][0], extracted_seqs[0], read[0][2], extracted_quals[0]]))) + '\n')
                elif pattern2:
                    # paired end sequencing, umi extracted from read2
                    newreads = [list(map(lambda x: x.strip(), [readnames[0], read[0][1], read[0][2], read[0][3]]))]
                    newreads.append(list(map(lambda x: x.strip(), [readnames[1], seqs[0], read[1][2], quals[0]])))
                    if keep_extracted and r2_extracted:
                        r2_extracted.write('\n'.join(list(map(lambda x: x.strip(), [read[1][0], extracted_seqs[0], read[1][2], extracted_quals[0]]))) + '\n')

                # write new reads to output fastq
                for i in range(len(outfastqs)):
                    outfastqs[i].write('\n'.join(newreads[i]) +'\n')
        else:
            NonMatching += 1
            if keep_discarded:
                for i in range(len(discarded_fastqs)):
                    discarded_fastqs[i].write('\n'.join(list(map(lambda x: x.strip(), read[i]))) + '\n')

    # close all open files
    for i in infastqs + outfastqs + discarded_fastqs + extracted_fastqs:
        i.close()
        
    print(Total, Matching, NonMatching)

    
if __name__ == '__main__':
        
    # create parser
    parser = argparse.ArgumentParser(prog='barcodex.py', description='A package for extracting Unique Molecular Identifiers (UMIs) from single or paired read sequencing data')
    subparsers = parser.add_subparsers()
       		
    # extract commands
    e_parser = subparsers.add_parser('extract', help="Extract UMIs from read sequences")
    e_parser.add_argument('--r1_in', dest='r1_in', help='Path to input FASTQ 1', required=True)
    e_parser.add_argument('--r1_out', dest='r1_out', help='Path to output FASTQ 1', required=True)
    e_parser.add_argument('--pattern', dest='pattern', help='Barcode string of regex for extracting UMIs in read 1')
    e_parser.add_argument('--pattern2', dest='pattern2', help='Barcode string of regex for extracting UMIs in read 2')
    e_parser.add_argument('--inline', dest='inline_umi', action='store_true', help='UMIs inline with reads or not. True if activated')
    e_parser.add_argument('--data', dest='data', choices=['single', 'paired'], default='single', help='Paired or single end sequencing')
    e_parser.add_argument('--r2_in', dest='r2_in', help='Path to input FASTQ 2. Fastq 2 for paired end sequencing with inline UMIs. Fastq with UMIs for single end sequencing with UMIs not in line')
    e_parser.add_argument('--r2_out', dest='r2_out', help='Path to output FASTQ 2')
    e_parser.add_argument('--r3_in', dest='r3_in', help='Path to input FASTQ 3. Fastq with UMIs for paired end sequencing with UMIs not in line')
    e_parser.add_argument('--separator', dest='separator', default='_', help='String separating the UMI sequence in the read name')
    e_parser.add_argument('--keep_extracted', dest='keep_extracted', action='store_true', help='Output the extracted UMIs and potentially discarded sequences from reads in separate fastqs. True if activated')
    e_parser.add_argument('--keep_discarded', dest='keep_discarded', action='store_true', help='Output reads with non-matching patterns to separate fastqs. True if activated')
    e_parser.add_argument('--full_match', dest='full_match', action='store_false', help='Requires the regex pattern to match the entire read sequence. True if activated')
    
    args = parser.parse_args()
    extract_barcodes(args.r1_in, args.r1_out, pattern=args.pattern, pattern2=args.pattern2,
                     inline_umi=args.inline_umi, data=args.data, keep_extracted=args.keep_extracted, keep_discarded=args.keep_discarded,
                     r2_in=args.r2_in, r2_out=args.r2_out, r3_in=args.r3_in, full_match=args.full_match)
        

    

    