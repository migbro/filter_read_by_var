#!/usr/bin/env python
"""
Adds variant info from variant report to results from filter_read_by_var.py
Usage: ./annotate_filter_results.py <report> <results>

Arguments:
<report>  Variant pipeline report file
<results> Result output with cross-species variant hits

Options:
-h

"""
from docopt import docopt
import sys

args = docopt(__doc__)

rpt = open(args['<report>'])
res = open(args['<results>'])

var_dict = {}
head = next(rpt).rstrip('\n')
# read through report
for line in rpt:
    info = line.rstrip('\n').split('\t')
    if len(info) > 1:
        if info[0] not in var_dict:
            var_dict[info[0]] = {}
        var_dict[info[0]][info[1]] = line.rstrip('\n')
rpt.close()
# print header
sys.stdout.write('Num in other species\tchrom in other species\tpos in other species\tread1 count\tread2 count'
                 '\tnum reads paired\t' + head + '\n')
for line in res:
    info = line.rstrip('\n').split('\t')
    try:
        sys.stdout.write('\t'.join(info[4:]) + '\t' + var_dict[info[0]][info[1]] + '\n')
    except:
        sys.stderr.write(line)
res.close()