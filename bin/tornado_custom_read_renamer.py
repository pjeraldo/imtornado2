#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 13:41:48 2013
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

parser= argparse.ArgumentParser(description="Renames reads according to a file with the read mappings.")

parser.add_argument('from_to_file', help="Tab-delimited file with from -> to mappings")
parser.add_argument('in_fasta', help="Input fasta file.")
parser.add_argument('out_fasta', help="Output fasta file with IDs renamed")

args= parser.parse_args()

trim= lambda s: s.strip().split()
id_map= {trim(l)[0]:trim(l)[1] for l in open(args.from_to_file, "rU")}

#renamed reads
renamed= (SeqRecord(seq=r.seq, id=id_map[r.id], description="original_id={}".format(r.id)) for r in SeqIO.parse(args.in_fasta, "fasta"))

#Save those reads... order doesn't matter if the mapping is complete
SeqIO.write(renamed, args.out_fasta, "fasta")