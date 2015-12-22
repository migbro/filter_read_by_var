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


def mutect_check(read_obj, align_file, skip_dict):
    frac = 0.3
    if read_obj.qname in skip_dict:
        return skip_dict, 0
    try:
        slen = 0
        clip = 0
        ins = []
        ins_size = []
        clip_pos = []
        clip_size = []
        m = re.findall('(\d+\w)', read_obj.cigarstring)
        for cig in m:
            slen += int(cig[:-1])
            # need to track number of bases clipped AND move backwards of clipped
            if cig[-1] == 'S' or cig[-1] == 'H':
                clip += int(cig[:-1])
                clip_pos.append(slen)
                clip_size.append(int(cig[:-1]))

            # need to move position later if there's an insertion forward
            if cig[-1] == 'I':
                ins.append(slen)
                ins_size.append(int(cig[:-1]))
        if not read_obj.is_proper_pair:
            sys.stderr.write('Not a proper pair!\n')

        if float(clip)/slen < frac and read_obj.mapping_quality > 0:
            if read_obj.has_tag('MD'):
                mm = re.findall('(\d+\w)', read_obj.get_tag('MD'))
                rqual = 0
                cur = 0
                for md in mm:
                    pos = cur + int(md[:-1])
                    for i in xrange(0, len(clip_pos), 1):
                        if pos > clip_pos[i]:
                            pos -= clip_size[i]
                            break
                    for i in xrange(0, len(ins), 1):
                        if pos > ins[i]:
                            pos += ins_size[i]
                            break
                    rqual += read_obj.query_alignment_qualities[pos]
                    cur = pos + 1
                if rqual > 100:
                    #pdb.set_trace()
                    sys.stderr.write('Sum quality of mismatches too high, skipping\n')
                    return skip_dict, 0
            if abs(read_obj.tlen) < 202:
                align_obj = pysam.AlignmentFile(align_file, 'rb')
                test = align_obj.mate(read_obj)
                if test.mapping_quality > read_obj.mapping_quality:
                    align_obj.close()
                    return skip_dict, 0
                else:
                    skip_dict[test.qname] = 1
                    align_obj.close()
                    return skip_dict, 1

            else:
                skip_dict[read_obj.qname] = 1
                return skip_dict, 1
    except:
        # pdb.set_trace()
        sys.stderr.write('Error!\n')
        return skip_dict, 0

hsa_file = pysam.AlignmentFile(args['<bam1>'])
bam_file = pysam.AlignmentFile(hsa_file, 'rb')
vcf_file = pysam.VariantFile(args['<vcf>'], 'r')
# output reads hitting variant in file
read_out = open('var_read_hits.txt', 'w')
i = 0
reads = {}
# store all variant objects
var_objs = []
# flag whether variant found in mouse
var_flag = {}
j = 1
mod = 100
err_ct = 0

to_skip = {}
for variant in vcf_file.fetch():
    if j % mod == 0:
        sys.stderr.write('Processing variant ' + str(j) + ' in vcf file\n')
    j += 1

    var_objs.append(variant)

    var = variant.alts[0]
    for read in bam_file.fetch(variant.chrom, variant.pos, (variant.pos + 1)):
        pos = variant.pos - read.pos - 1
        if pos >= 0 and read.is_proper_pair:
            try:
                m = re.findall('(\d+\w)', read.cigarstring)
                for cig in m:
                    if cig[-1] == 'D':
                        pos -= int(cig[:-1])
                read.query_alignment_sequence[pos]
            except:
                err_ct += 1
                continue
            # first check that read matches variant, then if it'd pass mutect filters
            if read.query_alignment_sequence[pos] == var:
                flag = 0
                (to_skip, flag) = mutect_check(read, hsa_file, to_skip)
                if flag == 1:
                    reads[read.qname] = {}
                    reads[read.qname]['pos'] = pos
                    read_out.write('\t'.join((read.qname, bam_file.getrname(read.tid), str(variant.pos))) + '\n')
                    reads[read.qname]['var'] = var
                    reads[read.qname]['v_idx'] = i
                    # pdb.set_trace()
    i += 1
sys.stderr.write(str(err_ct) + ' reads skipping in bam1 due to invalid cigar or positioning\n')
bam_file.close()
vcf_file.close()
read_out.close()

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

mmu_subset_bam = pysam.AlignmentFile(mmu_filtered, 'rb')
# sys.stderr.write('Indexing filtered reads\n')
# mmu_subset_bai = pysam.IndexedReads(mmu_subset_bam, 1)
# mmu_subset_bai.build
for read in mmu_subset_bam.fetch():
#    pdb.set_trace()
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
                var_flag[index]['chr'] = mmu_subset_bam.getrname(read.tid)
                var_flag[index]['pos'] = read.pos + cur_pos
                var_flag[index]['r1'] = 0
                var_flag[index]['r2'] = 0
                var_flag[index]['paired'] = 0
                var_flag[index]['var'] = 0
            var_flag[index]['var'] += 1
            if read.is_read1:
                var_flag[index]['r1'] += 1
            elif read.is_read2:
                var_flag[index]['r2'] += 1
            if read.is_paired:
                var_flag[index]['paired'] += 1


    except:
        # sys.stderr.write('Error at read ' + str(j) + ' skipping!\n')
        err_ct += 1
        continue
mmu_bam.close()
sys.stderr.write(str(err_ct) + ' mouse reads skipped due to missing cigar or invalid positioning\n')
out = open(args['<out>'], 'w')
for index in var_flag:
    out.write('\t'.join(
            (var_objs[index].chrom, str(var_objs[index].pos), var_objs[index].ref,
             var_objs[index].alts[0])) + '\t' + '\t'.join((str(var_flag[index]['var']), var_flag[index]['chr'],
                                                           str(var_flag[index]['pos']), str(var_flag[index]['r1']),
                                                           str(var_flag[index]['r2']),
                                                           str(var_flag[index]['paired']))) + '\n')
