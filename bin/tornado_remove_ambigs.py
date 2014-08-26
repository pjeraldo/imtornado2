#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 15:44:51 2013
"""
__author__ = "Patricio Jeraldo"
__copyright__ = "Copyright 2014, Mayo Foundation for Medical Education and Research "
__credits__ = ["Patricio Jeraldo", "Nicholas Chia"]
__license__ = "GPL"
__version__ = "2.0.0"
__maintainer__ = "Patricio Jeraldo"
__email__ = "jeraldo.patricio@mayo.edu"
__status__ = "Production"

import re
from Bio import SeqIO
import argparse

parser= argparse.ArgumentParser(description="Remove reads that have any ambiguous nucleotides (WSMKRYBDHVN).")

parser.add_argument('in_fasta', help="Input fasta file.")
parser.add_argument('out_fasta', help="Output fasta file with no ambigous nucleotides in reads.")

args= parser.parse_args()

#regular expression matching any ambiguous nucleotides,
#upper and lower case
ambigs= re.compile('[WSMKRYBDHVNwsmkrybdhvn]+')

#load sequences, retain if there are no ambigs in it
clean= (r for r in SeqIO.parse(args.in_fasta, "fasta") if ambigs.search(str(r.seq)) == None)

#save reads
SeqIO.write(clean, args.out_fasta, 'fasta')


