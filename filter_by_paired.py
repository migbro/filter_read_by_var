#!/usr/bin/env python
'''
Usage: ./filter_by_paired.py <bam1> <out>

Arguments:
<bam1> bam file to extract reads from
<out>  out file prefix

Options:
-h

'''
import re
import sys
import pysam
from docopt import docopt

args = docopt(__doc__)
hsa_file = args['<bam1>']
bam_file = pysam.AlignmentFile(hsa_file, 'rb')
out_filtered = args['<out>']
out = pysam.AlignmentFile(out_filtered, 'wb', template=bam_file)
pair_err = 0
cig_err = 0
for read in bam_file:
    if read.is_proper_pair:
        try:
            m = re.findall('(\d+\w)', read.cigarstring)
            out.write(read)
        except:
            cig_err += 1
    else:
        pair_err +=1
out.close()
bam_file.close()
pysam.index(out_filtered)
sys.stderr.write(str(pair_err) + ' unpaired reads, ' + str(cig_err) + ' missing cigar string\n')
