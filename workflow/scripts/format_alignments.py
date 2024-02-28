import argparse
import sys
import shutil
import os
import re


usage = """

"""

description = """
Simple python script that will remove putative contaminant sequences from an alignment based on strings in the sequence names. 
"""

# argparse
parser = argparse.ArgumentParser(usage=usage, description=description)
parser.add_argument("--input",     help = "Input directory containg '.fasta' files to be filtered.", nargs = "*", required=True)
parser.add_argument("--cont",      help = "Comma separated list of contaminant sample names.",       nargs = "*", required=False)
parser.add_argument("--output",    help = "Output directory to write output fasta files.",           required=True)
parser.add_argument("--overwrite", help = "Overwrite output directory.", action = "store_true",      required=False)
args = parser.parse_args()

# functions
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
        s_escape = re.escape(s)
        if re.search(s_escape, name):
            result = True
            print(f'Removing {name}')
            break
    return result

#assert str_present('Zet_ZKP_1315332_contig0;3595-5130;+;cox1', ['Zet_ZKP_1315332_contig0;3595-5130;+;cox1',  'Turbo_SRR15496837_contig0;12141-13676;-;cox1']) == True
#assert str_present('Zet_ZKP_1315332_contig0;3595-5130;+;cox1', ['Zet_ZKP_1315332',                           'Turbo_SRR15496837_contig0;12141-13676;-;cox1']) == True
#assert str_present('Zet_ZKP_1315332_contig0;3595-5130;+;cox1', ['cox1',                                      'Turbo_SRR15496837_contig0;12141-13676;-;cox1']) == True
#assert str_present('Zet_ZKP_1315332_contig0;3595-5130;+;cox1', ['Ilang_IWHI_ECP_1_contig0;3409-4944;+;cox1', 'Turbo_SRR15496837_contig0;12141-13676;-;cox1']) == False

# main 

print("\nRunning remove_contaminants.py")

# print input used
print("\nSearching for '.fasta' files in the following path(s):")
for i in args.input:
    print(f"   {i}")
#print(args.input)

if args.cont != None:
    print("\nRemoving sequences with names containing the following strings:")
    for i in args.cont:
        print(f"   {i}")
else:
    print("\nNo names provided to indicate putative contaminants")
#print(args.cont)

# create output dir
if os.path.exists(args.output):
    if args.overwrite: 
        print("\nOverwriting output directory")
        shutil.rmtree(args.output)
        os.mkdir(args.output)
    else: 
        sys.exit(f"Error: output directory {args.output} already exists.")
else: 
    print("\nCreating output directory")
    os.mkdir(args.output)

# iterate through fasta files in the input directory
for path in args.input:
    # print(path)
    for file in os.listdir(path):
        #print(file)
        if file.endswith(".fasta"):
            # read fasta
            print(f"\nReading fasta {file}")
            fas = read_fasta(f"{path}/{file}")
            # initiate tmp list to hold fasta data
            tmp = []
            # iterate throuh sequences in fasta
            for i in fas:
                name, seq = i[0], i[1]
                # write sequence to output list if contaminant not present in name
                if args.cont != None: 
                    if not str_present(name, args.cont):
                        name = format_name(name)
                        tmp.append([name, seq])
                else: 
                    name = format_name(name)
                    tmp.append([name, seq])

            # do not generate fasta if less than 5 sequences present
            if len(tmp) < 5: 
                print(f"No output generated for {file}. Need at least 5 sequences.")
            else:
                out = open(f"{args.output}/{file}", 'w')
                for x in tmp: 
                    name, seq = x[0], x[1]
                    out.write(f'>{name}\n{seq}\n')
                out.close()

print("Complete!")
