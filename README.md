# BarcodEX #

BarcodEX is tool for extracting Unique Molecular Identifiers (UMIs) from single or paired-end read sequences.
It can handle UMIs inline with reads or located in separate fastqs.

## Installation ##
### From PyPi ###
BarcodEx is available from PyPi

```pip install barcodex```

### From GitHub ###
Clone the BarcodEx repository

```git clone https://github.com/oicr-gsi/barcodex```

Install required python libraries by running
```pip install -r requirements.txt```

## Extraction of UMI sequences ##

BarcodEx extracts UMIs in read sequences or in fastqs with ```inline``` or ```separate``` sub-commands 

```barcodex -h
usage: barcodex.py [-h] [--umilist UMILIST] --prefix PREFIX
                   [--separator SEPARATOR]
                   {inline,separate} ...

A package for extracting Unique Molecular Identifiers (UMIs) from single or
paired read sequencing data

positional arguments:
  {inline,separate}     sub-command help
    inline              Extract UMIs located in read sequences
    separate            Extract UMIs located in separate file

optional arguments:
  -h, --help            show this help message and exit
  --umilist UMILIST     Path to file with valid UMIs (1st column)
  --prefix PREFIX       Specifies the start of the output files and stats json
                        files
  --separator SEPARATOR
                        String separating the UMI sequence in the read name```

### Extracts UMIs in read sequences ###

usage: ```barcodex --prefix PREFIX --separator SEPARATOR --umilist UMILIST inline --r1_in FASTQ1 --r2_in FATSQ2 --pattern1 PATTERN1 --pattern2 PATTERN2 --full_match``` 

Parameters

| argument | purpose | required/optional                                    |
| ------- | ------- | ------------------------------------------ |
| --r1_in | Path(s) to FASTQ(s) containing read 1   | required                |
| --r2_in | Path(s) to FASTQ(s) containing read 2   | optional              |
| --pattern1 | pattern or regex for extracting UMIs in read 1   | optional              |
| --pattern2 | pattern or regex for extracting UMIs in read 2   | optional              |
| --prefix | Specifies the start of the output files | required              |
| --separator | String separating the UMI sequence in the read name   | required              |
| --umilist | Path to file with valid UMIs | optional              |
| --full_match | Requires the regex pattern to match the entire read sequence | optional              |


### Extracts UMIs in fastqs ###

usage: ```barcodex --prefix PREFIX --separator SEPARATOR --umilist UMILIST separate --r1_in FASTQ1 --r2_in FATSQ2 --ru_in UMIs --pattern1 PATTERN1 --full_match``` 

Parameters

| argument | purpose | required/optional                                    |
| ------- | ------- | ------------------------------------------ |
| --r1_in | Path(s) to FASTQ(s) containing read 1   | required                |
| --r2_in | Path(s) to FASTQ(s) containing read 2   | optional              |
| --ru_in | Path(s) to FASTQ(s) containing UMIs     | required              |
| --pattern1 | pattern or regex for extracting UMIs in read 1   | optional              |
| --prefix | Specifies the start of the output files | required              |
| --separator | String separating the UMI sequence in the read name   | required              |
| --umilist | Path to file with valid UMIs | optional              |
| --full_match | Requires the regex pattern to match the entire read sequence | optional              |


BarcodEx extracts UMIs using either a pattern sequence or a regular expression and appends the concatenated UMIs to the read name preceded by a separator string specified in the command. 
UMIs can be extracted from read 1 and/or read 2 using respectively ```--pattern1``` and ```--pattern2```. At leat 1 pattern must be used. When extracting UMIs in read 1 and read 2, ```--pattern1``` and ```--pattern2``` must be both either a string sequence or a regular expression.
Reads that are not matching the provided patterns are discarded. Discarded reads are written to file for inspection. Morover, the extracted sequences are also recovered and written to file.
```--prefix``` specifies the start of the output files. (e.g ```--prefix``` x will result in x_R1.fastq.gz, x_R1.discarded.fastq.gz, x_R1.extracted.fastq.gz, x_UMI_counts.json, x_extraction_metrics.json) 


### Extraction with a string pattern ###

The pattern sequence must include one or more Ns, indicating the UMI bases, optionally followed by any nuccleotides corresponding to spacer sequence. 
For instance the pattern ```NNNNN``` extracts the first 5 nucleotides from the read whereas pattern ```NNNNNATCG``` extracts the first 9 nucleotides, appends nucleotides 1-5 to the read name and discard spacer ```ATCG```. Reads not matching ```NNNNNATCG``` are discarded. 

Extraction with the pattern sequence always extracts UMIs at the beginning of the read sequence. Extraction with regular expression offers more flexibility in the UMI design (see below).

As an example, consider read:

```
@MISEQ753:39:000000000-BDH2V:1:1101:17521:1593 1:N:0:
TCATGTCTGCTAATGGGAAAGAGTGTCCTAACTGTCCCAGATCGTTTTTTCTCACGTCTTTTCTCCTTTCACTTCTCTTTTTCTTTTTCTTTCTTCTTCTT
+
1>1A1DDF11DBDGFFA111111D1FEEG31AD1DAA1110BA00000//01A2A/B/B/212D2111D1222D12122B1B01D1@101112@D2D12BB
```

Extraction with pattern ```NNNNNNNNNNNNATGGGAAAGAGTGTCC``` will extract UMI ```TCATGTCTGCTA``` and add it to the read name. Spacer sequence ```ATGGGAAAGAGTGTCC``` is removed from read.
So the new read is now:

```
@MISEQ753:39:000000000-BDH2V:1:1101:17521:1593_TCATGTCTGCTA 1:N:0:
TAACTGTCCCAGATCGTTTTTTCTCACGTCTTTTCTCCTTTCACTTCTCTTTTTCTTTTTCTTTCTTCTTCTT
+
G31AD1DAA1110BA00000//01A2A/B/B/212D2111D1222D12122B1B01D1@101112@D2D12BB
```

Extracted sequence ```TCATGTCTGCTAATGGGAAAGAGTGTCC``` and its corresponding qualities ```1>1A1DDF11DBDGFFA111111D1FEE``` can be written to fastq file with option ```--keep_extracted``` (see below).

### Extraction with a regular expression ###

Regular expressions allow more flexibility for extracting UMIs, in particular UMIs with complex design and UMIs not starting at the beginning of the read.
A good introduction to regular expression can be found in this [Regular Expression HOWTO](https://docs.python.org/3/howto/regex.html). 
BarcodEx depends on the ```regex``` module rather than the standard ```re``` module because the former allows fuzzy matching.

Sequences are extracted from the read using named groups within the regex. Allowed named groups are ```umi``` and ```discard```. Syntax with named groups is as follow:
```(?<umi>.{3})(?<discard>T{2})```: extracts a 3bp UMI followed by TT spacer that is removed from read and discarded
The ```discard``` group removes nucleotides and qualities from the read while the ```umi``` group extracts the UMI that gets added to the read name.
Any sequence not contained in ```umi``` and ```discard``` groups will remain in the read. Thus, it is important to construct the regular expression such that the begining of the read is captured in groups.

For instance, consider the following read:

```
@MISEQ753:39:000000000-BDH2V:1:1101:17521:1593 1:N:0:
AATCGTCCATCG
+
1>1A1DDF11DB
```

The regex ```(?<umi>.{3})(?<discard>C{2})``` will extract UMI ```CGT``` and discard spacer ```CC```. But the first 3 nucleotides ```AAT``` will remain in the read with new read being:

```
@MISEQ753:39:000000000-BDH2V:1:1101:17521:1593_CGT 1:N:0:
AATATCG
+
1>111DB
```

To prevent the ```AAT``` before the UMI ```CGT``` to be part of the read, we need to account for the nucleotides upstream in the UMI in the regex ```(?<discard1>^.*)(?<umi>.{3})(?<discard2>C{2})```

```
@MISEQ753:39:000000000-BDH2V:1:1101:17521:1593_CGT 1:N:0:
ATCG
+
11DB
```

The extracted sequences and qualities can be recovered with option ```--keep_extract``` which writes the following read to file (see below).

```
@MISEQ753:39:000000000-BDH2V:1:1101:17521:1593_CGT 1:N:0:
AATCGTCC
+
1>1A1DDF
```

Multiple ```umi``` and ```discard``` named groups are allowed within the regex but they should be named differently. Naming is not important as long as groups contain the strings ```umi``` and ```discard```without special characters.
For instance the following 2 regex will give the same output:
- ```(?<discard_1>^.*)(?<umi>.{3})(?<discard_a>C{2})')``` and ```(?<discard1>^.*)(?<umi>.{3})(?<discard2>C{2})')```
- ```(?P<umi_1>^[ACGT]{3}[ACG])(?P<discard_1>T)|(?P<umi_2>^[ACGT]{3})(?P<discard_2>T)``` and ```(?P<umi_a>^[ACGT]{3}[ACG])(?P<discard_a>T)|(?P<umi_b>^[ACGT]{3})(?P<discard_b>T)```

### Filtering extracted UMIs against a list ###

The UMIs will only be accepted if they match an allow list provided with --umilist.
The list is a text file with one UMI per line. In the case of 2 reads with embedded UMIs, the two parts of the UMI must be on separate lines, optionally followed by the read number they apply to.
So, AAA would be allowed for either read 1 or read 2, while CCC 2 will allow CCC only on read 2.
It's also possible to write AAA 1 2 or AAA 1 and AAA 2 if desired.
 

### Extraction of UMIs in single or paired read sequences ### 

Parameters ```--r1_in``` and ```--r2_in``` indicate the paths to the read 1 and read 2 fastqs. ```--r1_in``` is always required and omitting ```--r2_in``` indicates single end sequences. 
BarcodEx requires gzipped input fastqs and outputs gzipped fastqs. Input fastqs for paired end data must be in sync.
```--prefix``` specifies the start of the output files (e.g. foo/bar/x_R1.fastq.gz, foo/bar/x_R2.fastq.gz)
 

### Extraction from multiple input fastqs ###

Multiple input fastqs can be processed together for read 1 and/or read 2 but generating a single output fastq for single end data and 2 output fastqs for paired read data.
The files mut be passed to ```--r1_in``` for read 1 fastqs and ```--r2_in``` for read 2 fastqs, each file being separated by white space.
The number of input fastqs for paired data must be the same for read 1 and read 2 and each list of files must be in the same order.

### Extraction of UMIs not inline with reads ###

Sub-command ```inline``` extracts UMIs located within the read sequences while sub-command ```separate``` is used to extract UMIs located in fastq.
With UMIs in file, ```--pattern1``` is used to extract UMIs from file(s) specified by ```--ru_in```. 

### Recovering discarded reads and extracted sequences ###

Reads without a matching pattern are written to file for inspection. ```--prefix``` specifies the start of the output files (e.g. foo/bar/x_R1.discarded.fastq.gz).
Extracted read sequences (UMIs and any spacer sequence removed from read, along with their qualities) are also written as a fastq file (e.g. foo/bar/x_R1/R2.extracted.fastq.gz for inline UMIs, and x_extracted.R2/3.fastq for UMIs located in files).
The fastqs with extracted sequences can be used with the main output fastqs to re-generate the original fastqs.


### UMIs and reads metrics ###

Two files with metrics of the UMI extraction are written in json format in the same directory specified by ```--prefix```.
```--prefix``` foo/bar/x results in  foo/bar/x_extraction_metrics.json and foo/bar/x_UMI_counts.json.
The first json captures information about the extraction process:
- total reads/pairs processed
- number of reads/pairs with matching pattern
- number of reads/pairs with non-matching pattern
- number of reads/pairs discarded due to unknown UMI
- pattern1
- pattern2
- umi-list

The second json records the UMI counts after extraction. For paired end data it counts the concatenated sequences with UMIs from read 1 and read 2, but adding a "." separator to track the read origin of each UMI.
For instance, "AAA.TCG": 10 in the json file indicates that sequence AAATCG is found 10 times in all the extracted UMIs, and that it is made of AAA from read 1 and TCG from read 2. One can then easily obtain counts of all UMIs from read 1 and read 2.


### Importing Barcodex as a module ###

BarcodEx can be run as a script or imported as a module to perform extraction within your own script.
The recommended import is ```from barcodex import extract_barcodes_inline``` or ```from barcodex import extract_barcodes_separate```. Dependent modules must also be imported (see Example command #8 below).


```
from barcodex import extract_barcodes_inline

help(barcodex.extract_barcodes_inline)
Help on function extract_barcodes_inline in module barcodex:

extract_barcodes_inline(r1_in, pattern1, prefix, pattern2=None, r2_in=None, full_match=False, separator='_', umilist=None)
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
```


### Example commands ###

**Example 1. Paired end end with string sequence**

Extraction of UMIs in read 1 for paired end data with inline UMIs with a string sequence.
Extracts the first 12 bp UMI when followed by spacer ATGGGAAAGAGTGTCC and remove spacer from read.
Umis are preceded by an underscore in the read header.

```
python3.6 barcodex.py \
--separator "_"
--prefix output \
inline
--r1_in myfile_R1.fastq \
--r2_in myfile_R2.fastq \
--pattern NNNNNNNNNNNNATGGGAAAGAGTGTCC \
```

**Example 2. Paired end with regex**

Extraction of UMIs in read 1 and read 2 for paired end data with inline UMIs with a regular expression.
Extracts the first 3 nucleotides as UMI and discards the next 2 nucleotide spacer sequence from each read.
Umis are preceded by an underscore in the read header.

```
python3.6 barcodex.py \
--separator "_"
--prefix output \
inline
--r1_in myfile_R1.fastq.gz \
--r2_in myfile_R2.fastq.gz \
--pattern1 "(?<umi>.{3})(?<discard>.{2})" \
--pattern2 "(?<umi>.{3})(?<discard>.{2})" \
```

**Example 3. Full match regex option**

Same as Example 2, but the regex patterns are modified to suit the ```--full_match``` regex requirement.

```
python3.6 barcodex.py
--separator "_"
--prefix output \
inline
--r1_in myfile_R1.fastq.gz \
--r2_in myfile_R2.fastq.gz \
--pattern1 "(?<umi>.{3})(?<discard>.{2}.+)" \
--pattern2 "(?<umi>.{3})(?<discard>.{2}.+)" \
--full_match
```

**Example 4. List of UMIs**

Extraction of UMIs in read 1 and read 2 for paired end data with inline UMIs with a regular expression.
Extracts a 4 bp UMI not ending with T and discard a following T spacer or a 3 bp UMI and discard the following T spacer.
UMIs start at the beginning of the read sequence and only higher caps A, T, C and G are allowed.
Extracted UMIs are checked against the true_barcode.txt UMI list. Reads with non-valid UMIs (ie. not present in true_barcode.txt) are discarded.
Umis are preceded by an underscore in the read header.

```
python3.6 barcodex.py \
--separator "_"
--prefix output \
--umilist true_barcodes.txt \
inline
--r1_in myfile_R1.fastq.gz \
--r2_in myfile_R2.fastq.gz \
--pattern1 "(?P<umi_1>^[ACGT]{3}[ACG])(?P<discard_1>T)|(?P<umi_2>^[ACGT]{3})(?P<discard_2>T)" \
--pattern2 "(?P<umi_1>^[ACGT]{3}[ACG])(?P<discard_1>T)|(?P<umi_2>^[ACGT]{3})(?P<discard_2>T)"
```

**example 5. Single end**

Extraction of UMIs in read 1 for single end with inline UMIs with a regular expression.
Extracts the first 12 bp UMI when followed by spacer ATGGGAAAGAGTGTCC and remove spacer from read.
Umis are preceded by an underscore in the read header.

```
python3.6 barcodex.py \
--separator "_"
--prefix x \
--umilist true_barcodes.txt \
inline
--r1_in myfile_R1.fastq.gz \
--pattern1 "(?P<umi>^.{12})(?P<discard>ATGGGAAAGAGTGTCC)" \
```

**Example 6. Multiple input fastqs**

Extraction of UMIs in read 1 and read 2 for paired end data with inline UMIs with a regular expression from 4 input fastq1 and 4 input fastq2.
Extracts the first 3 nucleotides as UMI and discards the next 2 nucleotide spacer sequences from each read.
Umis are preceded by an underscore in the read header.

```
python3.6 barcodex.py \
--separator "_"
--prefix output_umis \
--umilist true_barcodes.txt \
inline
--r1_in myfile_R1.1.fastq.gz myfile_R1.2.fastq.gz myfile_R1.3.fastq.gz myfile_R1.4.fastq.gz 
--r2_in myfile_R2.1.fastq.gz myfile_R2.2.fastq.gz myfile_R2.3.fastq.gz myfile_R2.4.fastq.gz 
--pattern1 "(?<umi>.{3})(?<discard>.{2})" \
--pattern2 "(?<umi>.{3})(?<discard>.{2})" \
```

**Example 7. Paired end with UMIs in file**

Extraction of UMIs in read 1 and read 2 for paired end data with UMIs located in read 3 with a regular expression.
Extracts the first 10 nucleotides as UMI. Umis are preceded by an underscore in the read header.

```
python3.6 barcodex.py \
--separator "_"
--prefix foo/bar/myfastq.umis \
separate
--r1_in myfile_R1.fastq.gz \
--r2_in myfile_R2.fastq.gz \ 
--ru_in myfile_R3.fastq.gz \ ## file with UMIs
--pattern1 "(?<umi>^.{10})" \
```

**Example 8. Importing BarcodEx as module within a script**

Sample as Example 2.

```
import gzip
from itertools import zip_longest
import regex
import json
import time
from barcodex import extract_barcodes_inline


extract_barcodes_inline(['myfile_R1.fastq.gz'], pattern1='(?<umi>.{3})(?<discard>.{2})', prefix='x',
                   pattern2='(?<umi>.{3})(?<discard>.{2})', r2_in=['myfile_R2.fastq.gz'], full_match=False, separator='_', umilist=None)

```
