#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 27 13:43:47 2014

"""

__author__ = "Patricio Jeraldo"
__copyright__ = "Copyright 2014-2016, Mayo Foundation for Medical Education and Research "
__credits__ = ["Patricio Jeraldo", "Nicholas Chia"]
__license__ = "GPL"
__version__ = "2.0.3.3"
__maintainer__ = "Patricio Jeraldo"
__email__ = "jeraldo.patricio@mayo.edu"
__status__ = "Beta"

import argparse
import sys

parser= argparse.ArgumentParser(description="Flattens out fastA files (one line per sequence)")
parser.add_argument('-i','--input_fasta', help="Input fasta file. Default= STDIN.", default=sys.stdin, type=argparse.FileType('r'))
parser.add_argument('-o','--output_fasta', help="Output fasta file. Default= STDOUT.", default=sys.stdout, type=argparse.FileType('w'))

args= parser.parse_args()

#Start with empty ID and sequence
current_header=''
current_seq=''

for line in args.input_fasta:
    my_line= line.strip()
    #if I'm the header line
    if my_line.startswith(">"):
        #If the seq is not empty:
        if current_header != '' and current_seq !='':
            #write out the current fasta record
            args.output_fasta.write(">" + current_header + "\n")
            args.output_fasta.write(current_seq + "\n")            
            #clear out the seq
            current_seq=''
        #Set the header as everything after the ">"
        current_header= my_line[1:]
    #else I'm part of the sequence
    else:
        current_seq+= my_line
#Finally, flush the last fasta record
args.output_fasta.write(">" + current_header + "\n")
args.output_fasta.write(current_seq + "\n")

#Done