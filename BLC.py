#!/usr/bin/env python3

import argparse
import re
from sys import argv

import ete3
from ete3 import Tree

from scipy import stats
from numpy import array
import numpy as np

import glob

###### Documentation

parser = argparse.ArgumentParser(description="Requires python 3 and the python package ete3. This script was written by Lily Hughes, lilychughes@gmail.com. ")
parser.add_argument('-g', '--genetrees' , dest = 'genetrees' , type = str , default= None , required= True, help = 'Directory with newick-formatted constrained ML gene trees.')
parser.add_argument('-g', '--ext' , dest = 'ext' , type = str , default= ".treefile" , required= False, help = 'File extension of newick-formatted constrained ML gene trees. Default .treefile')
parser.add_argument('-c', '--concat' , dest = 'concat' , type = str , default= None , required= True, help = 'Newick-formatted ML concatenated tree.')
parser.add_argument('-r', '--ratio' , dest = 'ratio' , type = int , default= 5 , required= False, help = 'Threshold ratio for branch lengths above which taxa are flagged. Default 5 (5X greater than concatenated branch length).')
parser.add_argument('-o', '--outgroups' , dest = 'outgroups' , type = str , default= None , required= True, help = 'Text file with the outgroup taxa. One taxon name per line.')

args, unknown = parser.parse_known_args()
#######

# import & root concatenated tree

c = Tree(args.concat)

og = open(args.outgroups)
cout = []
for line in og:
    line = line.strip("\n")
    cout.append(line)
og.close()
croot = c.get_common_ancestor(cout)
c.set_outgroup(croot)

# get a list of constrained gene trees to process
gtrees = []

for g in glob.glob(args.genetrees+"*"+args.ext):
    gtrees.append(g)


# FUNCTIONS    
def CommonTips(tree1, tree2):
'''
get a list of common tips between two trees
'''
    tree1Tips = tree1.get_leaf_names()
    tree2Tips = tree2.get_leaf_names()
    commonTips = []
    for tip in tree1Tips:
        if tip in tree2Tips:
            commonTips.append(tip)
    return(commonTips)
            
def GetBranchLengths(tree):
'''
this function returns a list of branch lengths (terminal and internal), minus the 0.0 branch length that ete3 puts in at the root
there should be n - 3 internal branch lengths + n terminal branch lengths where n is equal to the number of tips in the tree
'''
    blengths = []
    for node in tree.traverse("postorder"):
        parentNode = node.up
        blen = tree.get_distance(node, parentNode)
        blengths.append(blen)
    return blengths


def TerminalBranchLength(taxon, tree):
'''
this function gets the terminal branch length for a specified taxon and tree
'''
    parentNode = (tree&taxon).up
    terminalBL = tree.get_distance(taxon, parentNode)
    return(terminalBL)


for g in genetrees:
    t = Tree(g)

    common_taxa = CommonTips(t, c)
    
    t.prune(common_taxa)
    c.prune(common_taxa)
    
    
    # root the gene tree and order the nodes in the same way as in the concatenated tree
    c.ladderize()
    children = c.get_children()
    if children[0].is_leaf() == True:
        ancestor = children[0].get_leaf_names()[0]
        t.set_outgroup(ancestor)
    else:
        outgroups = children[0].get_leaf_names()
        ancestor = t.get_common_ancestor(outgroups)
        if t == ancestor:
            mp = t.get_midpoint_outgroup()
            t.set_outgroup(mp)        
            ancestor2 = t.get_common_ancestor(outgroups)
            t.set_outgroup(ancestor2)
        else:    
            t.set_outgroup(ancestor)
    t.ladderize()
      
    # make a list of terminals with ratios above the threshold
    ratioTaxa = []
    
    # calculate the ratios between trees for each taxon
    
    for taxon in common_taxa:
        ratio = TerminalBranchLength(taxon, t) / TerminalBranchLength(taxon, c)
        #print(taxon+":"str(ratio))
        if ratio >= 5:
            ratioTaxa.append(taxon)
    
    
    if len(ratioTaxa) > 0:
        out = open(g+".flagged_taxa.txt", "w")
        for name in ratioTaxa:
            out.write(name+"\n")
        out.close()    
    
    
    # prune the flagged taxa out of the tree
    
    if len(ratioTaxa) > 0:
        keepTaxa = []
        for leaf in c.get_leaf_names():
            if leaf not in ratioTaxa:
                keepTaxa.append(leaf)
        t.prune(keepTaxa)
        c.prune(keepTaxa)
    
    
    # make lists of branch lengths for each tree
    t_bl = GetBranchLengths(t)    
    c_bl = GetBranchLengths(c)
    
    
    # convert to array format for regression
    x = array(t_bl)
    y = array(c_bl)
    
    
    # Get R2 between the trees
    
    slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    
    # Write it to a file
    out = open(g+".R2.txt", "w")
    out.write("R-squared for "+g+": "+str(r_value**2))
    out.close()
    