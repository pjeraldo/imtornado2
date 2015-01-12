#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 11:27:40 2013
"""

__author__ = "Patricio Jeraldo"
__copyright__ = "Copyright 2014, Mayo Foundation for Medical Education and Research "
__credits__ = ["Patricio Jeraldo", "Nicholas Chia"]
__license__ = "GPL"
__version__ = "2.0.3"
__maintainer__ = "Patricio Jeraldo"
__email__ = "jeraldo.patricio@mayo.edu"
__status__ = "Production"

import argparse
from Bio import SeqIO

parser= argparse.ArgumentParser(description="Selects reads from a fasta file using a list of sequence IDs.")

parser.add_argument('id_list', help="File with sequence IDs to be selected")
parser.add_argument('in_fasta', help="Input fasta file.")
parser.add_argument('out_fasta', help="Output fasta file.")

args= parser.parse_args()

#Get list of IDs
accnos= [l.strip() for l in open(args.id_list, "rU")]

#a bit slower, but it guarantees a certain order
r_index= SeqIO.index(args.in_fasta, "fasta")

selected= (r_index[i] for i in accnos)

#Save those reads
SeqIO.write(selected, args.out_fasta, "fasta")