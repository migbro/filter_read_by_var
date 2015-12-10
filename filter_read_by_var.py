#!/usr/bin/env python
'''
Usage: ./filter_read_by_var.py <vcf> <bam1> <bam2> <out>

Arguments:
<vcf> vcf file with variant to search for
<bam1> bam file to extract reads from
<bam2> bam file to check reads in
<out>  out file prefix

Options:
-h

'''

import pysam
from docopt import docopt

args = docopt(__doc__)

bam_file = pysam.AlignmentFile(args['<bam1>'], 'rb')
vcf_file = pysam.VariantFile(args['<vcf>'], 'r')
i = 0
reads = []
var_objs = []
var_flag = {}
for variant in vcf_file.fetch():
    var_objs.append(variant)
    var_flag[i] = 0
    var = variant.alts[0]
    for read in bam_file.fetch(variant.chrom, variant.pos):
        pos = variant.pos - read.pos - 1
        if read.query_alignment_sequence[pos] == var:
            reads[read.qname] = []
            reads[read.qname]['pos'] = pos
            reads[read.qname]['var'] = var
            reads[read.qname]['v_idx'] = i
            # pdb.set_trace()
    i += 1
bam_file.close()
vcf_file.close()

mmu_bam = pysam.AlignmentFile(args['<bam2>'], 'rb')

for read in mmu_bam:
    if read.qname in reads and read.seq[reads[read.qname]['pos']] == reads[read.qname]['var']:
        var_flag[reads[read.qname]['v_idx']] += 1
mmu_bam.close()
out = open(args['<out>'], 'w')
for index in var_flag:
    if var_flag[index] > 0:
        out.write('\t'.join(
                (var_objs[index].chrom, str(var_objs[index].pos), var_objs[index].ref, var_objs[index].alts[0])) + str(
                var_flag[index]) + '\n')
