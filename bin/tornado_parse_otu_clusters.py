#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 13 15:33:29 2013
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


parser=argparse.ArgumentParser(description="Parses a usearch cluster to extract members of OTUs and failed to be assigned to any OTU.")

parser.add_argument('otu_fasta', help="Fasta file with OTU representatives." )
parser.add_argument('uc_file', help="USEARCH cluster file")
parser.add_argument('out_otus', help="Output file with OTUs.")
parser.add_argument('out_failures', help="Output with reads that failed to be assigned to an OTU.")

args= parser.parse_args()

#read the actual OTU representative IDs from the OTU file.
otu_ids= [int(r.id) for r in SeqIO.parse(args.otu_fasta, "fasta")]
#sort the ids
otu_ids.sort()
#create the dict to hold the OTU members
otus={rep:[] for rep in otu_ids}
#also create an empty list to hold the rejected IDs
failures=[]

#no go over the clusters file and parse it
for line in open(args.uc_file, "rU"):
    fields= line.strip().split('\t')
    if fields[0] == 'H':
        member_id= fields[8].split(' ')[0]
        #I assume it's numeric
        rep_id= int(fields[9].split(' ')[0])
        #add member to otus
        otus[rep_id].append(member_id)
    if fields[0] == 'N':
        member_id= fields[8].split(' ')[0]
        #add rejected to failures list
        failures.append(member_id)

#save the OTUs
with open(args.out_otus, "w") as outfile:
    for otu in otus.iterkeys():
        outfile.write("{}\t".format(otu))
        outfile.write("\t".join(otus[otu]) + '\n')

#save the failures
#sort them
failures.sort()
with open(args.out_failures, "w") as outfile:
    outfile.write('\n'.join(failures) + '\n')

#Done!        


