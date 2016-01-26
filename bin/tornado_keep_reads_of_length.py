#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
"""
__author__ = "Patricio Jeraldo"
__copyright__ = "Copyright 2014, Mayo Foundation for Medical Education and Research "
__credits__ = ["Patricio Jeraldo", "Nicholas Chia"]
__license__ = "GPL"
__version__ = "2.0.3"
__maintainer__ = "Patricio Jeraldo"
__email__ = "jeraldo.patricio@mayo.edu"
__status__ = "Alpha"

import argparse
from Bio import SeqIO

parser= argparse.ArgumentParser(description="Discards reads without the required length.")

parser.add_argument('trim_length', help="Read length required to keep the read. Anything less and the reads is gone.", type=int)
parser.add_argument('in_fasta', help="Input fasta file.")
parser.add_argument('out_fasta', help="Output file with non-conforming reads removed.")

args= parser.parse_args()

good= (r for r in SeqIO.parse(args.in_fasta, "fasta") if len(r) == args.trim_length)

SeqIO.write(good, args.out_fasta, "fasta")
