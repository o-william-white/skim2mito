
import argparse
import sys
from ete3 import Tree
import re

# argparse
parser = argparse.ArgumentParser()
parser.add_argument("--input",    help = "Input newick tree.",                 required=True)
parser.add_argument("--output",   help = "Output newick file.",                required=True)
parser.add_argument("--outgroup", help = "Outgroup(s) e.g. 'sampleX,sampleY'", required=True)
args = parser.parse_args()

# read input tree
tree = Tree(open(args.input, "r").read().rstrip("\n"))

# print info
print(f"Input tree:\n{tree}\n")
print(f"Rooting on:\n{args.outgroup}\n")

# split outgroups
outgroups_tmp = args.outgroup.split(",")

# function to get leaf name if outgroup is not an exact match
# returns None if no match found
def get_outgroup_leaf_name(tree, outgroup):
    for n in tree.traverse():
        if n.is_leaf():
            if re.search(outgroup, n.name):
                return n.name

# iterate through outgroups to check if matches are exact or partial
outgroups = []
for o in outgroups_tmp: 
    if o == get_outgroup_leaf_name(tree, o):
        print(f"Exact match found for {o}") 
        outgroups.append(o)
    else: 
        partial_match = get_outgroup_leaf_name(tree, o)
        if partial_match != None:
            print(f"Partial match found for {o}: {partial_match}")
            outgroups.append(partial_match)
        else: 
            print(f"No match found for {o}")

# function to get first ingroup sample
def select_first_ingroup_leaf(tree, outgroups):
    leaf_names = []
    for n in tree.traverse():
        if n.is_leaf() and n.name not in outgroups:
            leaf_names.append(n.name)
    return(leaf_names[0])

# root tree

# if single outgroup providied, find node
if len(outgroups) == 1:
    node = tree&outgroups[0]
    print("Setting root")
    tree.set_outgroup(node)
    tree.ladderize()
    print(tree)
else:
    # find ancestor of multiple samples
    print(f"HERE: {outgroups}")
    ancestor = tree.get_common_ancestor(outgroups)
    # check if ancestor is already the current root
    if tree == ancestor:
        print("MRCA of outgroups is the current root node")
        # check if current root is monophyletic
        is_mono = tree.check_monophyly(values = outgroups, target_attr="name")[0]
        if is_mono is True:
            print("Outgroups monophyletic. Will try to ladderise tree")
            tree.ladderize()
            print(tree)
        # root non-monophyletic but mrca is the the current root - paraphyletic
        else:
            print("Outgroups non-monophyletic. Will try rooting on another node before re-rooting")
            # root on another node
            tree.set_outgroup(select_first_ingroup_leaf(tree, outgroups))
            # root on outgroup ancestor again
            ancestor = tree.get_common_ancestor(outgroups)
            tree.set_outgroup(ancestor)
            tree.ladderize()
            print(tree)
    else:
        print("Setting root")
        tree.set_outgroup(ancestor)
        tree.ladderize()
        print(tree)

# write outfile
print()
print("Writing output")
outfile = open(args.output, "w")
outfile.write(tree.write())
outfile.close()

