#!/usr/bin/env python
'''
Usage: ./mutect_param_bam_filter.py <bam1>

Arguments:
<bam1> bam file to extract reads from


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
out = pysam.AlignmentFile('clip_filtered.bam', 'wb', template=bam_file)
frac = 0.3
to_skip = {}

for read in bam_file.fetch():
    if read.qname in to_skip:
        continue
    try:
        slen = 0
        clip = 0
        dflag = 0
        ins = []
        ins_size = []
        clip_pos = []
        clip_size = []
        m = re.findall('(\d+\w)', read.cigarstring)
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
        if not read.is_proper_pair:
            sys.stderr.write('Not a proper pair!\n')

        if float(clip)/slen < frac and read.mapping_quality > 0 and read.is_proper_pair:
            if read.has_tag('MD'):
                mm = re.findall('(\d+\w)', read.get_tag('MD'))
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
                    rqual += read.query_alignment_qualities[pos]
                    cur = pos + 1
                if rqual > 100:
                    #pdb.set_trace()
                    sys.stderr.write('Sum quality of mismatches too high, skipping\n')
            if abs(read.tlen) < 202:
                test = bam_file.mate(read)
                if test.mapping_quality > read.mapping_quality:
                    out.write(test)
                else:
                    out.write(read)
                to_skip[test.qname] = 1
            else:
                out.write(read)
                to_skip[read.qname] = 1
    except:
        # pdb.set_trace()
        sys.stderr.write('Cigar missing!\n')
bam_file.close()
out.close()
pysam.sort('clip_filtered.bam', 'filtered_sorted')
pysam.index('filtered_sorted.bam')
