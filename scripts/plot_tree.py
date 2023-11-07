import argparse
from ete3 import Tree, NodeStyle, TreeStyle, TextFace

import os
os.environ["DISPLAY"] = "8888"

# from xvfbwrapper import Xvfb

# argparse
parser = argparse.ArgumentParser()
parser.add_argument("--input",  help = "Input tree in newick format.",    required=True)
parser.add_argument("--output", help = "Output plot in png format.",      required=True)
parser.add_argument("--size",   help = "Output plot size in mm (width).", required=True)
args = parser.parse_args()

# read tree file
tree_file = open(args.input, "r")
newick = tree_file.readline()
tree_file.close()

# create tree object for ete3
t = Tree(newick)

# tree style
ts = TreeStyle()
ts.show_leaf_name = False
ts.min_leaf_separation = 20
ts.margin_top = 20
ts.show_scale = False
ts.draw_guiding_lines = True

## function to add tip label
def format_labels(node):
    # if node is tip
    if node.is_leaf():
        # add text face
        label_face  = TextFace(node.name, fsize=10, fgcolor="black")
        node.add_face(label_face,  column=0, position="aligned")

# apply function to all nodes
for n in t.traverse():
    format_labels(n)

# function to add bootstrap values as text on top of the branches
def add_bootstrap_values(node):
    # if node has a bootstrap value and is not root or tip
    if node.support and not node.is_leaf() and not node.is_root():
        # add the bootstrap value as text face on top of the branch
        bootstrap_face = TextFace(int(node.support), fsize=10, fgcolor="darkgray")
        node.add_face(bootstrap_face, column=0, position="branch-top")
        
# apply function to all nodes
for n in t.traverse():
    add_bootstrap_values(n)

# function to format nodes
def format_node(node):
    # if node is a leaf node
    if node.is_leaf():
        leaf_node = NodeStyle()
        leaf_node["size"] = 2
        leaf_node["fgcolor"] = "black"
        node.set_style(leaf_node)
    else:
        # normal node
        internal_node = NodeStyle()
        internal_node["size"] = 0
        node.set_style(internal_node)
    node.hz_line_width = 100
    node.vt_line_width = 100
        
# apply function to all nodes
for n in t.traverse():
    format_node(n)

#vdisplay = Xvfb()
#vdisplay.start()

t.render(args.output, w=args.size, units="mm", tree_style=ts)

#vdisplay.stop()

