#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 13 16:43:31 2013
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
import sys
import biom
#check version of the BIOM API
biom_version= biom.__version__
table_from_API={}
if biom_version[0] == '1':
    biom_API= '1.0'
elif biom_version[0] == '2':
    biom_API= '2.1'
else:
	sys.stderr.write("Unsupported BIOM-format API: {}\n".format(biom_version))
	sys.exit(1)
if biom_API == '2.1':
    from biom.table import Table
    table_from_API[biom_API]=Table
else:
    from biom.table import table_factory
    table_from_API[biom_API]=table_factory
import numpy
from itertools import izip


parser= argparse.ArgumentParser(description="Creates a BIOM file from the otus found using usearch/uparse.")

parser.add_argument('otu_file', help="File with OTUs.")
parser.add_argument('mapping_file', help="File with sample metadata. Must be QIIME-conformant.")
parser.add_argument('taxonomy_file', help="File with the taxonomy calls associated with an OTU. Can be consensus or taxonomy of the representative read")
parser.add_argument('out_biom_file', help="Output BIOM file with OTU table and metadata")

args= parser.parse_args()

#split by tab
def st(s): return s.strip().split('\t')

#open the mapping file
with open(args.mapping_file, 'rU') as infile:
    headers= st(infile.next())
    
    samples= []
    metadata= []
    for line in infile:
        if line.strip() is '' or line[0] == '#':
            continue
        data= st(line)
        assert len(data) == len(headers)
        samples.append(data[0])
        meta= {h:d for h,d in izip(headers[1:], data[1:])}
        metadata.append(meta)

sample_dict= {name:i for i,name in enumerate(samples)}
#make taxonomy dict
taxonomy_dict= {st(line)[0]:st(line)[1].split(';')[:-1] for line in open(args.taxonomy_file)}

#make otu dict
otu_dict= {st(line)[0]:st(line)[1:] for line in open(args.otu_file)}
    
otu_id_dict= {id:i for i,id in enumerate(otu_dict.iterkeys())}
#and make the appropriate metadata in the correct order
taxonomy_metadata= [{"taxonomy":taxonomy_dict[otu]} for otu in otu_dict.iterkeys()]

#make the observations matrix
# number_of_otus x number_of_samples
observations= numpy.zeros((len(otu_dict.keys()),len(samples)), dtype=int)

#now fill the otu thing
for otu,reads in otu_dict.iteritems():
    for read in reads:
        sample_name= read.split('_')[0]
        sample_id= sample_dict[sample_name]
        otu_id= otu_id_dict[otu]
        #add one to each count
        observations[otu_id,sample_id]+= 1
#now that the observation matrix is filled, create the biom table

otu_table= table_from_API[biom_API](data=observations, sample_ids= samples, observation_ids= otu_dict.keys(), sample_metadata= metadata, observation_metadata= taxonomy_metadata)

#Write the biom file
with open(args.out_biom_file, "w") as outfile: 
    if biom_API=='2.1':
        outfile.write(otu_table.to_json("IM-TORNADO-{}".format(__version__)) + '\n')
    else:
        outfile.write(otu_table.getBiomFormatJsonString("IM-TORNADO-{}".format(__version__)) + '\n')

#done
