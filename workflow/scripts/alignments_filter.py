
import argparse

# argparse
parser = argparse.ArgumentParser()
parser.add_argument("--input",     help = "Input fasta",  required=True)
parser.add_argument("--output",    help = "Output fasta", required=True)
parser.add_argument("--threshold", help = "Maximum percentage of missing data", required=True, type=float)
args = parser.parse_args()

# read fasta
def read_fasta(filename):
    name, seq = None,""
    fasta = open(filename)
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

# count missing
def filter_seq(fasta, threshold, output):
    output_fasta = open(output, "w")
    for i in fasta:
        name, seq = i[0], i[1]
        seq_length = len(seq)
        seq_missing = seq.count("-") 
        pct_missing = seq_missing / seq_length
        print(f"{name}\t{seq_length}\t{seq_missing}\t{pct_missing}")
        if pct_missing < threshold:
            output_fasta.write(f">{name}\n{seq}\n")
    output_fasta.close()
    
# read fasta
fas = read_fasta(args.input)

# filter fasta
filter_seq(fas, args.threshold, args.output)

