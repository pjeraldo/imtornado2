#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 12 16:43:23 2012
"""
__author__ = "Patricio Jeraldo"
__copyright__ = "Copyright 2014-2016, Mayo Foundation for Medical Education and Research "
__credits__ = ["Patricio Jeraldo", "Nicholas Chia"]
__license__ = "GPL"
__version__ = "2.0.3.3"
__maintainer__ = "Patricio Jeraldo"
__email__ = "jeraldo.patricio@mayo.edu"
__status__ = "Production"

import argparse

parser= argparse.ArgumentParser(description="Finds the read IDs common to two FASTA formatted file.")

parser.add_argument("fasta_1", help="First FASTA file.", type=argparse.FileType('rU'))
parser.add_argument("fasta_2", help="Second FASTA file.", type=argparse.FileType('rU'))
parser.add_argument("output", help="Output file with the read IDs common to both FASTA files.", type=argparse.FileType('w'))

args= parser.parse_args()

names_1= {line.strip().strip(">") for line in args.fasta_1 if line.strip()[0] is '>' and line.strip()[0] is not ''}
names_2= {line.strip().strip(">") for line in args.fasta_2 if line.strip()[0] is '>' and line.strip()[0] is not ''}

common= {name for name in names_1 if name in names_2}

args.output.write('\n'.join(common) + '\n')
args.output.close()
   