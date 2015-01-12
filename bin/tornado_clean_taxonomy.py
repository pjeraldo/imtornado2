#!/usr/bin/env python
"""
Created on Tue Jan 21 15:36:59 2014
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
import argparse

parser= argparse.ArgumentParser(description="Filter unneeded lines from taxonomy file")


parser.add_argument('in_fasta', help="Input fasta file with OTUs to be kept.")
parser.add_argument('in_taxa', help="Taxonomy file with entries to be removed.")
parser.add_argument('out_taxa', help="Output taxonomy file.")

args= parser.parse_args()

#Read IDs from fasta
otu_ids= {r.id for r in SeqIO.parse(args.in_fasta, "fasta")}

#Lines to be kept
clean_taxa= (line.strip() for line in open(args.in_taxa, "rU") if line.strip().split()[0] in otu_ids)

#Write everything
with open(args.out_taxa, "w") as outfile:
    outfile.write('\n'.join(clean_taxa) + '\n')
    
#Done!