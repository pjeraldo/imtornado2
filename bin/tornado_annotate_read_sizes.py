#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 10:54:41 2013
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
from Bio.SeqRecord import SeqRecord

parser= argparse.ArgumentParser(description="Renames IDs in the OTU file.")

parser.add_argument('cutoff', help="Length cutoff for the reads.", type=int)
parser.add_argument('in_fasta', help="Input fasta file with unique sequences.")
parser.add_argument('in_counts', help="Input counts file with number of replicates of the representative. See mothur's count.seqs")
parser.add_argument('out_fasta', help="Output file with size annotations in the format expeted by usearch.")

args= parser.parse_args()

#get the counts file
prep= lambda s: s.strip().split()

counts= {prep(l)[0]:prep(l)[1] for l in open(args.in_counts, "rU")}

annotated= (SeqRecord(seq=r.seq, id="{0};size={1};".format(r.id, counts[r.id]), description=r.description) for r in SeqIO.parse(args.in_fasta, "fasta") if len(r) >= args.cutoff)

SeqIO.write(annotated, args.out_fasta, "fasta")
