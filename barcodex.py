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

    :param filename (str): File name or file path    
    
    Return True if the file is gzipped and False otherwise
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
    
    :param fastq (str): Path to the file fastq (compressed or not)
    
    Return an opened fastq file with file handler at the begining of the file
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
    
    :param readname: read header
    :param UMI: UMI sequence
    :param separator: string separating the UMI sequence and part of the read header
    
    Returns the read name with UMI sequence
    '''
    
    readname = readname.split(' ')
    readname[0] = readname[0] + separator + UMI
    readname = ' '.join(readname)
    return readname


def _is_pattern_sequence(pattern):
    '''
    (str) -> bool
    
    :param pattern: Pattern to be extracted from reads
    
    Return True if all elements of pattern are valid nucleotides
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
    (str, str) -> None
    
    '''
    
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
    

def _extract_from_sequence(read, UMI, spacer, inline_umi=True):
    '''
    (list, str, str, bool) -> (str, str, str, str, str)
    
    Returns a tuple with the read sequence and qualities after barcode extraction,
    the umi sequence, the extracted read sequence and qualities
    '''
    
    # initialize variables
    seq, qual, umi_seq, extracted_seq, extracted_qual = '', '', '', '', ''
    
    # extract UMI starting at begining of the read
    if spacer in read[1]:
        if spacer == read[1][len(UMI): len(UMI) + len(spacer)]:
            umi_seq = read[1][: len(UMI)]
            # no extraction if umi not inline with read
            if inline_umi:
                extracted_seq = read[1][: len(UMI) + len(spacer)]
                extracted_qual = read[3][: len(UMI) + len(spacer)]
                seq = read[1][len(UMI) + len(spacer):]
                qual = read[3][len(UMI) + len(spacer):]
    
    return seq, qual, umi_seq, extracted_seq, extracted_qual
    


def _extract_from_regex(read, p, full_match=False, inline_umi=True):
    '''
    (list, _regex.Pattern, bool, bool)
    
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
        umi_seq = ''
        for i in umi_pos:
            umi_seq += read[1][i[0]:i[1]]
        # no extraction if umi not inline with read
        if inline_umi:
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



def _check_fastq_sync(read1, read2, read3 = None):
    '''
    (str, str, str | None) -> None
    
    
    '''
    
    L = [read1, read2]
    if read3:
        L.append(read3)
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
                raise ValueError('Expecting paired end sequences with inline UMIs. \
                                 Paths to r1 and r2 I/O fastqs required')
            if not r3 is None:
                raise ValueError('Expecting paired end sequences with inline UMIs. \
                                 Paths to r1 and r2 I/O fastqs required. Path to r3 input not needed')
        elif not inline_umi:
            # requires r1, r2 and r3
            # requires pattern, pattern2 not needed
            if any(map(lambda x: x is None, [r1_in, r1_out, r2_in, r2_out, r3_in])):
                raise ValueError('Expecting paired end sequences with out of read UMIs. \
                                 Paths to r1 and r2 I/O fastqs and to r3 input fastq required')
    elif data == 'single':
        if inline_umi:
            # requires r1, but not r2 and not r3
            if any(map(lambda x: x is None, [r1_in, r1_out])):
                raise ValueError('Expecting single end sequences with inline UMIs. \
                                 Paths to r1 I/O fastqs required')
            if any(map(lambda x: x is not None, [r2_in, r2_out, r3_in])):
                raise ValueError('Expecting single end sequences with inline UMIs. \
                                 Paths to r1 I/O fastqs required. Paths to r2 I/O fastqs and r3 input fastq not needed')
        elif not inline_umi:
            # requires r1 and r2, r3 not needed
            if any(map(lambda x: x is None, [r1_in, r1_out, r2_in])):
                raise ValueError('Expecting single end sequences with out of read UMIs. \
                                 Paths to r1 I/O and r2 input fastqs required')
            if any(map(lambda x: x is not None, [r2_out, r3_in])):
                raise ValueError('Expecting single end sequences with out of read UMIs. \
                                 Paths to r2 output and r3 input fastq not needed')


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
        
    
    
    
    
    

def extract_barcodes(r1_in, r1_out, pattern, pattern2=None, inline_umi=True,
                     data='single', keep_extracted=True, keep_discarded=True,
                     r2_in=None, r2_out=None, r3_in=None, full_match=False):
    """





    """
    
    # to do:
    # remove extension before appending discarded or extracted extension
    # check that gzip.open opens for compression
    # option to output fastq or fastq.gz
    # module
    # script
    # write stats file
    # umis in 5' or/and 3'
    # keep discarded and extracted in sync
    
    
    # check input and output parameters
    _check_input_output(r1_in, r1_out, data, inline_umi, r2_in, r2_out, r3_in)
    # check pattern parameters 
    _check_pattern_options(pattern, pattern2, data, inline_umi)

    # open files for reading
    r1 = _open_fastq(r1_in)    
    r2 = _open_fastq(r2_in) if r2_in else None
    r3 = _open_fastq(r3_in) if r3_in else None
    
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
        
        
    # check if pattern is nucleotide string or regex
    if IsPatternSequence(pattern) == True:
        seq_extract = True
        # get UMi and spacer
        UMI, spacer = get_umi_spacer(pattern)
    else:
        seq_extract = False
        # compile pattern
        p = regex.compile(pattern)
             
    # check if second pattern is used
    if pattern2:
        # check that both patterns are either strings or regex
        CheckExtractionMode(pattern, pattern2)
        if IsPatternSequence(pattern2) == True:
            # get UMi and spacer
            UMI2, spacer2 = get_umi_spacer(pattern2)
        else:
            # compile pattern
            p2 = regex.compile(pattern2)
    
         
    # do a check based on number of input and output files currently supported in the library prep ini
    if num_reads == 3:
        assert actual_reads == 2, 'Expecting 2 output fastqs and 3 input fastqs'
    elif num_reads == 2:
        assert actual_reads == 2, 'Expecting 2 output fastqs and 2 input fastqs'
    elif num_reads == 1:
        assert actual_reads == 1, 'Expecting 1 output fastq and 1 input fastq'

    print("Preprocessing reads...")
    
    # make a list of fastqs open for reading
    fastqs = [i for i in [r1, r2, r3] if i != None]
    
    # make a list of files open for writing
    writers = [i for i in [r1_writer, r2_writer] if i != None]

    # make a list of files open for writing discarded reads
    if spacer == True:
        discarded = [i for i in [r1_discarded, r2_discarded, r3_discarded] if i != None]
    
    # check the number of input files
    # create iterator with reads from each file
    if len(fastqs) == 3:
        I =  zip(getread(fastqs[0]), getread(fastqs[1]), getread(fastqs[2]))
    elif len(fastqs) == 2:
        I =  zip(getread(fastqs[0]), getread(fastqs[1]))
    elif len(fastqs) == 1:
        I =  zip(getread(fastqs[0]))
    
    # count all reads and reads with incorrect and correct umi/spacer configuration
    Total, Correct, Incorrect = 0, 0, 0
    # Record all umi sequences with correct umi/spacer configuration    
    UmiSequences = []
    # Record read length
    ReadLength = [set(), set()]
        
    # loop over iterator with slices of 4 read lines from each line
    for reads in I:
        # count total reads
        Total += 1
        # extract umi sequences from reads
        if seq_extract == True:
            # extraction using string sequence. assumes UMi starts at begining of read
            
            
            seq, qual, umi_seq, extracted_seq, extracted_qual = ExtractFromSequence(read, UMI, spacer)
        
        
        
        
            print('hello')
        
        # make a list of read sequences
        readseqs = [i[1] for i in reads]
        umis = extract_umis(readseqs, umi_locs, umi_lens, umi_pos)
        
        # make a list of reads with umis
        reads_with_umis = [readseqs[i-1] for i in umi_locs]
        
        # skip reads with spacer in wrong position
        if spacer == True and False in [correct_spacer(reads_with_umis[i], umis[i], umi_pos[i] -1, spacer_seq) for i in range(len(reads_with_umis))]:
            # count reads with incorrect umi/spacer configuration 
            Incorrect += 1
            # write reads with incorrect umi/spacer configuration to separate fastqs 
            for i in range(len(discarded)):
                for j in range(len(reads[i])):
                    discarded[i].write(reads[i][j])
            # skip read processing
            continue

        # count number of reads with correct umi/spacer configuration
        Correct += 1
        
        # collect umi sequences
        UmiSequences.extend(umis)
        
        # edit read names and add umi
        # make parallel lists with begining and end of read name from r1, and from r2 or r3
        read_name1, rest1 = reads[0][0].rstrip().split(' ')
        readnames, namerests = [read_name1], [rest1]
        
        if len(reads) > 1:
            # edit read name from r2 (or r3)
            read_name2, rest2 = reads[-1][0].rstrip().split(' ')
            readnames.append(read_name2)
            namerests.append(rest2)
         
        # make lists with umi lengths, spacer lengths and umi positions for     
        UmiLength, SpacerLength, UmiPositions = [umi_len_r1, umi_len_r2], [spacer_len_r1, spacer_len_r2], [umi_pos_r1, umi_pos_r2]    
        
        for i in range(len(writers)):
            # if paired reads and umis are in each read: concatenate umis and assign concatenated umi to each read
            # if paired read and single umi: assign the same umi to each read
            # if single end read, assign the umi to its read
            if len(umis) > 1:
                # concatenate umis
                umiseq = ''.join(umis)
            elif len(umis) == 1:
                umiseq = umis[0]
            # add umi to read name and write to outputfile
            writers[i].write(readnames[i] + ":" + umiseq + " " + namerests[i] + "\n")
            # remove umi and spacer from read seq. write remaining of read to outputfile
            if i == 0:
                # determine index for reads k <- 0 for r1, -1 for r2 or r3
                k = 0
            elif i > 0:
                k = -1
            
            # check if umis are inline with reads or not
            # remove umi + spacer (if any) if umis inline with reads
            if umi_inline == False:
                p, l, s = 0, 0, 0
            else:
                p, l, s = UmiPositions[i], UmiLength[i], SpacerLength[i]
            
            # compute read length
            ReadLength[i].add(len(reads[k][1][p + l + s:]))
            
            # write new fastqs
            writers[i].write(reads[k][1][p + l + s:])
            writers[i].write(reads[k][2])
            writers[i].write(reads[k][3][p + l + s:])
        
    # close all open files
    for i in writers:
        i.close()
    for i in fastqs:
        i.close()
    if spacer == True:
        if len(discarded) != 0:
            for i in discarded:
                i.close()
    
    print("Complete. Output written to {0}".format(outdir))
    
    # record read indo as json
    D = {'Total': Total, 'Correct': Correct, 'Incorrect': Incorrect}
    if len(ReadLength[0]) != 0:
        lr1 = list(ReadLength[0])[0]
        D['length_read1'] = lr1   
    if len(ReadLength[1]) != 0:
        lr2 = list(ReadLength[1])[0]
        D['length_read2'] = lr2
    return D, UmiSequences



if __name__ == '__main__':
        
    # create parser
    parser = argparse.ArgumentParser(prog='barcodex.py', description='A package for extracting Unique Molecular Identifiers (UMIs) from single or paired read sequencing data')
    subparsers = parser.add_subparsers()
       		
    # extract commands
    e_parser = subparsers.add_parser('extract', help="Preprocess mode for processing fastq files")
    p_parser.add_argument('-o', '--OutDir', dest='outdir', help='Output directory. Available from command or config')
    p_parser.add_argument('-r1', '--Read1', dest='read1', help='Path to first FASTQ file.', required=True)
    p_parser.add_argument('-r2', '--Read2', dest='read2', help='Path to second FASTQ file, if applicable')
    p_parser.add_argument('-r3', '--Read3', dest='read3', help='Path to third FASTQ file, if applicable')
    p_parser.add_argument('-p', '--Prepname', dest='prepname', choices=['HALOPLEX', 'SURESELECT', 'EPIC-DS', 'SIMSENSEQ-PE', 'SIMSENSEQ-SE'], 
                          help='Name of library prep to  use (defined in library_prep_types.ini)', required=True)
    p_parser.add_argument('-pf', '--Prepfile', dest='prepfile', help='Path to your library_prep_types.ini file')
    p_parser.add_argument('-c', '--Config', dest='config', help='Path to your config file')
    p_parser.add_argument('-px', '--Prefix', dest= 'prefix', help='Prefix for naming umi-reheradered fastqs. Use Prefix from Read1 if not provided') 
    #p_parser.set_defaults(func=extract_barcodes)
    



    
    
        
    args = parser.parse_args()
    try:
        args.func(args)
    except AttributeError as e:
        print(e)
        print(parser.format_help())
        
        
    extract_barcodes(r1_in, r1_out, pattern, pattern2=None, inline_umi=True,
                     data='single', keep_extracted=True, keep_discarded=True,
                     r2_in=None, r2_out=None, r3_in=None, full_match=False):

    
    ices=[True, False], type=ConvertArgToBool, default=False, help='Keep the most abundant family and ignore families at other positions within each group. Default is False')
        