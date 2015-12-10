#!/usr/bin/env python
'''
Usage: ./filter_read_by_var.py <vcf> <bam>

Arguments:
<vcf> vcf file with variant to search for
<bam> bam file to extract reads from

Options:
-h

'''
import sys
import pysam
import pdb
from docopt import docopt
args = docopt(__doc__)

bam_file = pysam.AlignmentFile(args['<bam>'], 'rb')
vcf_file = pysam.VariantFile(args['<vcf>'], 'r')

for variant in vcf_file.fetch():
    for read in bam_file.fetch(variant.chrom, (variant.pos - 100), (variant.pos + 100)):
        pdb.set_trace()
