#!/usr/bin/env python

__author__ = "Patricio Jeraldo"
__copyright__ = "Copyright 2014, Mayo Foundation for Medical Education and Research "
__credits__ = ["Patricio Jeraldo", "Nicholas Chia"]
__license__ = "GPL"
__version__ = "2.0.2"
__maintainer__ = "Patricio Jeraldo"
__email__ = "jeraldo.patricio@mayo.edu"
__status__ = "Production"

from Bio import SeqIO
import argparse

parser= argparse.ArgumentParser(description="Converts a FASTQ formatted file (qualities in PHRED33 format) into a FASTA formatted file, and optionally a file with corresponding qualities.")

parser.add_argument('fastq', help="Input FASTQ file.")
parser.add_argument('fasta', help="Output FASTA file.")
parser.add_argument('--qual', help="Optional QUAL file to store read qualities" )

args= parser.parse_args()

SeqIO.convert(args.fastq, "fastq", args.fasta, "fasta")

if args.qual is not None:
    SeqIO.convert(args.fastq, "fastq", args.qual, "qual") 
