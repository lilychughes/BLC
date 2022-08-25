#!/usr/bin/env python3

import argparse
import re
from sys import argv

import ete3
from ete3 import Tree

from scipy import stats
from numpy import array
import numpy as np

###### Documentation

parser = argparse.ArgumentParser(description="Requires python 3 and the python package ete3. This script was written by Lily Hughes, lilychughes@gmail.com. ")
parser.add_argument('-t', '--tree' , dest = 'tree' , type = str , default= None , required= True, help = 'Newick-formatted constrained ML gene tree.')
parser.add_argument('-c', '--concat' , dest = 'concat' , type = str , default= None , required= True, help = 'ROOTED newick-formatted ML concatenated tree.')
parser.add_argument('-o', '--output' , dest = 'output' , type = str , default= None , required= True, help = 'Name of output file to write taxa with terminal branch lengths above the threshold ratio.')
parser.add_argument('-r', '--ratio' , dest = 'ratio' , type = int , default= 5 , required= False, help = 'Threshold ratio for branch lengths above which taxa are flagged. Default 5.')

args, unknown = parser.parse_known_args()
#######

# some stuff to import trees and prune taxa here

t = Tree(args.tree)
c = Tree(args.concat)


# get a list of common tips between two trees
def CommonTips(tree1, tree2):
    tree1Tips = tree1.get_leaf_names()
    tree2Tips = tree2.get_leaf_names()
    commonTips = []
    for tip in tree1Tips:
        if tip in tree2Tips:
            commonTips.append(tip)
    return(commonTips)
            
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


# this function gets the terminal branch length for a specified taxon and tree
def TerminalBranchLength(taxon, tree):
    parentNode = (tree&taxon).up
    terminalBL = tree.get_distance(taxon, parentNode)
    return(terminalBL)

    

# make a list of terminals with ratios above the threshold
ratioTaxa = []

# calculate the ratios between trees for each taxon (I used the absolute value so that either tree can be the denominator)

for taxon in common_taxa:
    ratio = TerminalBranchLength(taxon, t) / TerminalBranchLength(taxon, c)
    #print(taxon+":"str(ratio))
    if ratio >= 5:
        ratioTaxa.append(taxon)


if len(ratioTaxa) > 0:
    out = open(args.output, "w")
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

# this function returns a list of branch lengths (terminal and internal), minus the 0.0 branch length that ete3 puts in at the root
# there should be n - 3 internal branch lengths + n terminal branch lengths where n is equal to the number of tips in the tree
def GetBranchLengths(tree):
    blengths = []
    for node in tree.traverse("postorder"):
        parentNode = node.up
        blen = tree.get_distance(node, parentNode)
        blengths.append(blen)
    return blengths


# make lists of branch lengths for each tree
t_bl = GetBranchLengths(t)    
c_bl = GetBranchLengths(c)


# convert to array format for regression
x = array(t_bl)
y = array(c_bl)


# Get R2 between the trees

slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)

print("R-squared: "+str(r_value**2))

