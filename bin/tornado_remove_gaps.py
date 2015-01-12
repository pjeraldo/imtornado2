#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue May 15 11:39:00 2012
"""
__author__ = "Patricio Jeraldo"
__copyright__ = "Copyright 2014, Mayo Foundation for Medical Education and Research "
__credits__ = ["Patricio Jeraldo", "Nicholas Chia"]
__license__ = "GPL"
__version__ = "2.0.3"
__maintainer__ = "Patricio Jeraldo"
__email__ = "jeraldo.patricio@mayo.edu"
__status__ = "Production"

from Bio import SeqIO
from Bio.Alphabet import generic_nucleotide
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from bitarray import bitarray
import string
import argparse

def remove_gapped_cols(seq_records):
  #create alignment
  length= len(str(seq_records[0].seq))
  #using bitarrays.... generate translation tables
  trans_table= string.maketrans(string.uppercase + string.lowercase + '.-', ''.join(['0' for i in xrange(52)]) + '11')
  gap_bitarray= length*bitarray('1')
  zero= bitarray('0')
  #translate gaps to '1', characters to '0'.. bitwise AND with
  #previous colums
  for record in seq_records:
      gap_bitarray &= bitarray(str(record.seq).translate(trans_table))

  #list of non-gap columns, the ones we want to keep
  nongaps= gap_bitarray.search(zero)
  for rec in seq_records:
    seq_s= str(rec.seq)
    new_s= ''.join([seq_s[i] for i in nongaps])
    yield SeqRecord(seq=Seq(new_s, generic_nucleotide), id=rec.id, name='', description='')  
    
#Get input files

parser= argparse.ArgumentParser(description="Removes gap-only columns from FASTA-formatted multiple sequence alignment files.")

parser.add_argument('input', help="Input FASTA file.", type=argparse.FileType('r'))
parser.add_argument('output', help="Output FASTA file.", type=argparse.FileType('w'))

args= parser.parse_args()

records= list(SeqIO.parse(args.input, 'fasta'))

ungapped_recs= remove_gapped_cols(records)

SeqIO.write(ungapped_recs, args.output, 'fasta')
