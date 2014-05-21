#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sat May 10 14:52:30 2014

@author: patricio
"""

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse

parser= argparse.ArgumentParser(description="Converts fastq to fasta, renaming read IDs according to the sample name.")

parser.add_argument('in_fastq', help="Input fastq file.")
parser.add_argument('out_fasta', help="Output fasta file.")

args=parser.parse_args()

#Sanity check
assert args.in_fastq != args.out_fasta

#Ger the sample name. It is everything before the underscore
sample=args.in_fastq.split('_')[0]

#Now load up the file, rename file, add old id as annotation
renamed =(SeqRecord(seq=r.seq, id="{0}_{1}".format(sample,n), description="old_id={}".format(r.id)) for n,r in enumerate(SeqIO.parse(args.in_fastq, 'fastq'),1))

#Output the renamed files

SeqIO.write(renamed, args.out_fasta, 'fasta')