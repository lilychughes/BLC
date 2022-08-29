#!/usr/bin/env python3

import argparse
import re
from sys import argv

import ete3
from ete3 import Tree

import glob

###### Documentation

parser = argparse.ArgumentParser(description="Requires python 3 and the python package ete3. This script was written by Lily Hughes, lilychughes@gmail.com. ")
parser.add_argument('-a', '--alignments' , dest = 'alignments' , type = str , default= None , required= True, help = 'Path to directory with FASTA gene alignments')
parser.add_argument('-c', '--constraint' , dest = 'constraint' , type = str , default= None , required= True, help = 'Tree file to prune to create constraint tree.')
parser.add_argument('-e', '--ext', dest = 'ext', type = str, default=".fasta", required= False, help = 'File extension of FASTA formatted alignment files. Default is .fasta')
args, unknown = parser.parse_known_args()
######


# Open ML tree
c = Tree(args.constraint, format=9)

def CommonTips(tree, taxlist):
    """Function to get a list of common tips between two trees"""
    tree1Tips = tree.get_leaf_names()
    commonTips = []
    for tip in tree1Tips:
        if tip in taxlist:
            commonTips.append(tip)
    return(commonTips)

# Make a list of alignments to process
alignments = []

for a in glob.glob(args.alignments+"*"+args.ext):
    alignments.append(a)

def GetAliTaxa(ali):
    """Makes a list of taxon names in a FASTA alignment"""
    tips = []
    for line in ali:
        if line.startswith(">"):
            line = line.strip(">")
            taxon = line.strip("\n")
            tips.append(taxon)
    return tips        


for a in alignments:
    ali = open(a)
    atips = GetAliTaxa(ali)
    ali.close()
    commonTaxa = CommonTips(c, atips)
    constraint = c
    constraint.prune(commonTaxa)
    constraint.write(format=9, outfile=a+".constraint")
    
