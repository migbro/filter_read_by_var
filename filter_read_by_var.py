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
import re
import sys
import pdb
import pysam
from docopt import docopt

args = docopt(__doc__)

bam_file = pysam.AlignmentFile(args['<bam1>'], 'rb')
vcf_file = pysam.VariantFile(args['<vcf>'], 'r')
i = 0
reads = {}
# store all variant objects
var_objs = []
# flag whether variant found in mouse
var_flag = {}
j = 1
mod = 100
err_ct = 0
for variant in vcf_file.fetch():
    if j % mod == 0:
        sys.stderr.write('Processing variant ' + str(j) + ' in vcf file\n')
    j += 1

    var_objs.append(variant)

    var = variant.alts[0]
    #    for read in bam_file.fetch(variant.chrom, (variant.pos - 100), (variant.pos + 100)):
    for read in bam_file.fetch(variant.chrom, variant.pos, (variant.pos + 1)):
        pos = variant.pos - read.pos - 1
        if pos >= 0:
            try:
                m = re.findall('(\d+\w)', read.cigarstring)
                for cig in m:
                    if cig[-1] == 'D':
                        pos -= int(cig[:-1])
                read.query_alignment_sequence[pos]
            except:
                err_ct += 1
                continue
            if read.query_alignment_sequence[pos] == var:
                reads[read.qname] = {}
                reads[read.qname]['pos'] = pos
                reads[read.qname]['var'] = var
                reads[read.qname]['v_idx'] = i
                # pdb.set_trace()
    i += 1
sys.stderr.write(str(err_ct) + ' reads skipping in bam1 due to invalid cigar or positioning\n')
bam_file.close()
vcf_file.close()

mmu_bam = pysam.AlignmentFile(args['<bam2>'], 'rb')
j = 1
mod = 10000000
err_ct = 0
# two passes - one to print a bam with just reads, another index and then process by read
mmu_filtered = 'relevant_reads.bam'
out = pysam.AlignmentFile(mmu_filtered, 'wb', template=mmu_bam)
for read in mmu_bam:
    if j % mod == 0:
        sys.stderr.write('At read ' + str(j) + ' in bam 2 file\n')
    j += 1
    if read.qname in reads:
        out.write(read)
out.close()
sys.stderr.write('Creating filtered reads bam index\n')
pysam.index(mmu_filtered)
# pdb.set_trace()
mmu_bam.close()
mmu_subset_bam = pysam.AlignmentFile(mmu_filtered, 'rb')
# sys.stderr.write('Indexing filtered reads\n')
# mmu_subset_bai = pysam.IndexedReads(mmu_subset_bam, 1)
# mmu_subset_bai.build
for read in mmu_subset_bam.fetch():
    pdb.set_trace()
    try:
        # make same adjustment above for deletion

        cur_pos = reads[read.qname]['pos']
        m = re.findall('(\d+\w)', read.cigarstring)
        for cig in m:
            if cig[-1] == 'D':
                cur_pos -= int(cig[:-1])
        if read.qname in reads and read.seq[cur_pos] == reads[read.qname]['var']:
            index = reads[read.qname]['v_idx']
            if index not in var_flag:
                var_flag[index] = {}
                # store corresponding mouse info
                var_flag[index]['chr'] = mmu_bam.getrname(read.tid)
                var_flag[index]['pos'] = read.pos + cur_pos
                var_flag[index]['r1'] = 0
                var_flag[index]['r2'] = 0
                var_flag[index]['paired'] = 0
                var_flag[index]['var'] = 0
            var_flag[index]['var'] += 1
            if read.is_read1:
                var_flag[index]['r1'] += 1
            else:
                var_flag[index]['r2'] += 1
            if read.is_paired:
                var_flag[index]['paired'] += 1


    except:
        # sys.stderr.write('Error at read ' + str(j) + ' skipping!\n')
        err_ct += 1
        continue

sys.stderr.write(str(err_ct) + ' mouse reads skipped due to missing cigar or invalid positioning\n')
out = open(args['<out>'], 'w')
for index in var_flag:
    out.write('\t'.join(
            (var_objs[index].chrom, str(var_objs[index].pos), var_objs[index].ref,
             var_objs[index].alts[0])) + '\t' + '\t'.join((str(var_flag[index]['var']), var_flag[index]['chr'],
                                                           str(var_flag[index]['pos']), str(var_flag[index]['r1']),
                                                           str(var_flag[index]['r2']),
                                                           str(var_flag[index]['paired']))) + '\n')
