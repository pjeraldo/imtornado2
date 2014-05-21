#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 14:52:13 2013
"""
__author__ = "Patricio Jeraldo"
__copyright__ = "Copyright 2014, Mayo Foundation for Medical Education and Research "
__credits__ = ["Patricio Jeraldo", "Nicholas Chia"]
__license__ = "GPL"
__version__ = "2.0.0"
__maintainer__ = "Patricio Jeraldo"
__email__ = "jeraldo.patricio@mayo.edu"
__status__ = "Production"

import argparse
from Bio import SeqIO

parser= argparse.ArgumentParser(description="Renames IDs in the OTU file.")

parser.add_argument('id_list', help="List of IDs to remove from fasta file.")
parser.add_argument('in_fasta', help="Input fasta file.")
parser.add_argument('out_fasta', help="Output file with selected IDs removed.")

args= parser.parse_args()

accnos= {l.strip() for l in open(args.id_list)}

good= (r for r in SeqIO.parse(args.in_fasta, "fasta") if r.id not in accnos)

SeqIO.write(good, args.out_fasta, "fasta")
