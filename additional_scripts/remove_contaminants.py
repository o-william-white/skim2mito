import argparse
import re

usage = """
TBC
"""

description = """
 
"""

# argparse
parser = argparse.ArgumentParser(usage=usage, description=description)
parser.add_argument("--input",  help = "Input fasta file.", required=True)
parser.add_argument("--cont",   help = "Comma separated list of contaminant sample names.", required=True)
parser.add_argument("--output", help = "Output fasta file.", required=True)
args = parser.parse_args()

# read fasta
def read_fasta(filename):
    name, seq = None,''
    fasta = open(filename, 'r')
    for line in fasta:
        if line.startswith('>') and name == None:
            name = line.rstrip('\n').replace('>','')
        else:
            if line.startswith('>') and name != None:
                yield [name, seq]
                name = line.rstrip('\n').replace('>','')
                seq = ''
            else:
                seq = seq + line.rstrip('\n')
    yield [name, seq]
    fasta.close()

def format_name(name):
    name = name.split(';')[0]
    name = re.sub("_contig\d*$|","", name)
    name = re.sub("_circular|","", name)
    return name

assert format_name('Spec_SPHI_MRT_3_contig0;6106-6801;+;atp6') == 'Spec_SPHI_MRT_3'
assert format_name('Spec_SPHI_MRT_3_circular;6106-6801;+;atp6') == 'Spec_SPHI_MRT_3'

def str_present(name, list_strings):
    result = False
    print(f'Checking {name}')
    for s in list_strings:
        if re.search(s, name):
            result = True
            print(f'Removing {name}')
    return result

#assert is_cont('Spec_SPHI_MRT_3', ['Suavo_61502', 'Zet_ZALP_103187', 'Suavo_61096']) == False
#assert is_cont('Spec_SPHI_MRT_3', ['Suavo_61502', 'Zet_ZALP_103187', 'Spec_SPHI_MRT_3']) == True

# read fasta
fas = read_fasta(args.input)

# samples to remove from fasta
samples = args.cont.split(',')

# output fasta
out = open(args.output, 'w')

for f in fas:
    name, seq = f[0], f[1]
    if not str_present(name, samples):
        name = format_name(name)
        out.write(f'>{name}\n{seq}\n')

# close output fasta
out.close()

