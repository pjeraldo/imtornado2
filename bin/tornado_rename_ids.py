#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 12:24:45 2013
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
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_nucleotide
from Bio.Seq import Seq
import argparse
import sys

def trim_reads_R1(l, reads, new_ids):  
    trim= lambda s: s[:l]
    for r in reads:
        if len(r) > l:
            yield SeqRecord(seq=Seq(trim(str(r.seq)),generic_nucleotide), id=new_ids[r.id], description="original_id={}".format(r.id))
        else:
            yield SeqRecord(seq=r.seq, id=new_ids[r.id], description="original_id={}".format(r.id))

def trim_reads_R2(l, reads, new_ids):   
    trim= lambda s: s[-l:]
    for r in reads:
        rc_seq= r.seq.reverse_complement()
        if len(r) > l:
            yield SeqRecord(seq=Seq(trim(str(rc_seq)),generic_nucleotide), id=new_ids[r.id], description="original_id={}".format(r.id))
        else:
            yield SeqRecord(seq=rc_seq, id=new_ids[r.id], description="original_id={}".format(r.id))


parser= argparse.ArgumentParser(description="Renames sequence ids to sample name. Reverse-complement R2. Output will be named PREFIX_R?.fasta. Also a new.groups file will be created for downstream use.")
parser.add_argument('--trim_r1', help='Maximum length of R1 reads. Default=250', type=int, default=250)
parser.add_argument('--trim_r2', help='Maximum length of R1 reads. Default=200', type=int, default=200)
parser.add_argument('PREFIX', help="PREFIX for the output files.")
parser.add_argument('in_fwd_reads', help="Input forward read files. Input ALL of them. It'll figure out the R2s. This is necessary due to numbering convetions.", nargs='*')

args= parser.parse_args()

#counter for id names
id_count=1
#out file names
out_fwd_file="{}_R1.fasta".format(args.PREFIX)
out_rev_file="{}_R2.fasta".format(args.PREFIX)
#for each fwd read
for fwd_file in args.in_fwd_reads:
    sample,suffix= fwd_file.split('_')
    new_suffix= suffix.replace('R1','R2')
    rev_file= "{0}_{1}".format(sample,new_suffix)
    sys.stderr.write("Processing {} and {}\n".format(fwd_file, rev_file))
    #get the current ids.. unfortunately I have to read all the files.
    #This is faster, in my experience (why? we're using more i/o this way)
    #Maybe there's more overhead in creating lists/sets
    r1_ids= {r.id for r in SeqIO.parse(fwd_file, "fasta")}
    r2_ids= {r.id for r in SeqIO.parse(rev_file, "fasta")}
    all_ids= r1_ids | r2_ids
    common_old_ids= r1_ids & r2_ids
    #construct new ids here to maintain consistency
    new_ids= {old:"{0}_{1}".format(sample,i) for i,old in enumerate(all_ids, id_count)}   
    common_new_ids= {new_ids[i] for i in common_old_ids}
    assert len(common_new_ids) <= len(new_ids)
    #no need to sort
    #process R1 first
    reads= SeqIO.parse(fwd_file, "fasta")
    new_reads= trim_reads_R1(args.trim_r1, reads, new_ids)
    #Save R1.
    SeqIO.write(new_reads, open(out_fwd_file, "a"), "fasta")
    #process R2 now
    reads= SeqIO.parse(rev_file, "fasta")
    new_reads= trim_reads_R2(args.trim_r2, reads, new_ids)
    #Save R2
    SeqIO.write(new_reads, open(out_rev_file, "a"), "fasta")
    #Append to a groups file
    with open("{}.groups".format(args.PREFIX), "a") as outfile:
        new=("{0}_{1}\t{0}".format(sample,i) for i in range(id_count,id_count+len(all_ids)))
        outfile.write('\n'.join(new) +'\n')
    #Write common IDs to a file
    with open("{}.common.accnos".format(args.PREFIX), "a") as outfile:
        outfile.write('\n'.join(common_new_ids) + '\n')

    id_count+= len(all_ids)

#save the last id in a file for use by orphan renamer

sys.stderr.write("Files {} and {} saved,\n".format(out_fwd_file, out_rev_file))
sys.exit(0)
    