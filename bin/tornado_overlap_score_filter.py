#!/usr/bin/env python

import numpy as n
from itertools import izip
import sys
import argparse

parser= argparse.ArgumentParser(description="Removes sequence reads that have poor alignment scores according to various criteria.")

parser.add_argument('input_scores', help="Input scores file generated when aligning reads using infernal")
parser.add_argument('output_bad_ids', help="Output files with list of read IDs with poor alignment scores.", type=argparse.FileType('w'))
parser.add_argument('-t','--tolerance', help='Tolerance for starting or ending position of a read. A read will be removed if its starting or endingposition differs in more than TOLERANCE bases with respect to the median starting or ending position, or if the alignment bitscore is negative. Default TOLERANCE value is 3.', type=int, default=3)

args= parser.parse_args()

in_fn= args.input_scores

accnos= [l.split()[1] for l in open(in_fn, 'r') if l[0] is not '#']
cmfrom= n.array([int(l.split()[3]) for l in open(in_fn, 'r') if l[0] is not '#'])
cmto= n.array([int(l.split()[4]) for l in open(in_fn, 'r') if l[0] is not '#'])
cmscore= n.array([float(l.split()[6]) for l in open(in_fn, 'r') if l[0] is not '#'])


start= n.median(cmfrom)
end= n.median(cmto)

sys.stderr.write("Median start position: {}\n".format(start))
sys.stderr.write("Median end position: {}\n".format(end))


bad_accnos= {acc for acc,cmf,cmt,score in izip(accnos,cmfrom,cmto,cmscore) if (abs(cmf - start) > args.tolerance) or (abs(cmt - end) > args.tolerance) or (score < 0)}

print("Total reads to remove: {}".format(len(bad_accnos)))

args.output_bad_ids.write('\n'.join(bad_accnos) + '\n')
args.output_bad_ids.close()
