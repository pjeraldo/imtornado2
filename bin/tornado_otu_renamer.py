#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 14:31:36 2013
"""
__author__ = "Patricio Jeraldo"
__copyright__ = "Copyright 2014, Mayo Foundation for Medical Education and Research "
__credits__ = ["Patricio Jeraldo", "Nicholas Chia"]
__license__ = "GPL"
__version__ = "2.0.2"
__maintainer__ = "Patricio Jeraldo"
__email__ = "jeraldo.patricio@mayo.edu"
__status__ = "Production"

import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

parser= argparse.ArgumentParser(description="Renames IDs in the OTU file.")

parser.add_argument('in_fasta', help="Input fasta file with otus.")
parser.add_argument('out_fasta', help="Output file with OTU IDs renamed.")

args= parser.parse_args()

renamed= (SeqRecord(seq=r.seq, id="{}".format(i), description="original_id={}".format(r.id)) for i,r in enumerate(SeqIO.parse(args.in_fasta, "fasta"),1))

SeqIO.write(renamed, args.out_fasta, "fasta")
