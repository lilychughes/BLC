#!/usr/bin/env python3

import re
from sys import argv
import argparse

import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import glob

###### Documentation

parser = argparse.ArgumentParser(description="Requires python3 and Biopython. Written by Lily Hughes lilychughes@gmail.com")
parser.add_argument('-a', '--alignments' , dest = 'alignments' , type = str , default= None , required= True, help = 'Directory with FASTA alignments to clean. This directory must ALSO contain the .flagged_taxa.txt files from BLC.py')
parser.add_argument('-e', '--ext', dest = 'ext', type = str, default=".fasta", required= False, help = 'File extension of FASTA formatted alignment files. Default is .fasta')
args, unknown = parser.parse_known_args()

#######

# Make a list of alignments to process
alignments = []

for a in glob.glob(args.alignments+"*"+args.ext):
    alignments.append(a)

# Prune taxa from alignments if they have been flagged by BLC.py
for ali in alignments:
    flag = glob.glob(ali+"*.flagged_taxa.txt")
    if len(flag) > 0:
        reference = open(ali)
        seq_dict = SeqIO.to_dict(SeqIO.parse(reference, "fasta"))
        reference.close()
        
        allTaxa = seq_dict.keys()
        
        toRemove = open(flag[0])
        
        TaxatoRemove = []
        
        for line in toRemove:
            taxon = line.strip("\n")
            TaxatoRemove.append(taxon)
        
        keepers =[]
        
        for item in allTaxa:
            if item not in TaxatoRemove:
                keepers.append(seq_dict[item])
        
        
        outfile = open(ali+".BLC_cleaned", "w")
        
        for entry in keepers:
            outfile.write(entry.format("fasta"))
        
        outfile.close()
        