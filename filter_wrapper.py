#!/usr/bin/env python
'''
Usage: ./filter_read_by_var.py <vcf_list> <bam1> <bam2>

Arguments:
<vcf_list> list of vcf files with variant to search for
<bam1> bam file to extract reads from
<bam2> bam file to check reads in

Options:
-h

'''
import pdb
import sys
import os
from docopt import docopt
sys.path.append('/home/ubuntu/TOOLS/Scripts/utility')
from job_manager import job_manager
th = 8
filter_py = '/home/ubuntu/TOOLS/filter_by_var/filter_read_by_var.py'

args = docopt(__doc__)
bam1 = args['<bam1>']
bam2 = args['<bam2>']
vcf_list = open(args['<vcf_list>'],'r')
cmd_list = []
for vcf in vcf_list:
    vcf = vcf.rstrip('\n')
    parts = os.path.basename(vcf).split('\.')
    out = parts[0] + '.hits.txt'
    cmd = ' '.join((filter_py, vcf, bam1, bam2, out))
    log_file = parts[0] + '.log'
    cmd += ' 2> ' + log_file
    cmd_list.append(cmd)
vcf_list.close()
job_manager(cmd_list, th)