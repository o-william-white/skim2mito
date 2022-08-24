import argparse
import requests
import time
import pandas as pd

# argparse

argparse_description = '''
Python script to annotate chloroplast assemblies using Chloe https://chloe.plastid.org/
'''

argparse_usage = '''
python chloe_api.py --input <input_fasta> --output <output_directory>
'''

parser = argparse.ArgumentParser(description=argparse_description, usage=argparse_usage)
parser.add_argument("--input", help="Input fasta to annotate", required=True)
parser.add_argument('--output', help="Output directory to write files to", required=True)
args = parser.parse_args()

# chloe
print('Submitting assembly')
fd = open(args.input, 'rb')
CHLOE = 'https://chloe.plastid.org'
res = requests.post(CHLOE + '/annotate', files=dict(file=fd)).json()

# allow time for job to complete
print('Running job')
time.sleep(60)

# get sff and gff outputs
sff = pd.read_csv(CHLOE + res['download_sff'], sep='\t', skiprows=1, header=None)
gff = pd.read_csv(CHLOE + res['download_gff'], sep="\t", skiprows=1, header=None,
                  names=["sequence_id", "source", "feature_type", "feature_start", "feature_end", "score", "strand",
                         "phase", "attributes"])

# format bed file
bed = gff[gff["feature_type"] == "gene"]
bed["gene"] = bed["attributes"].str.split(";").str[1].str.replace("Name=", "")
bed = bed[["sequence_id", "feature_start", "feature_end", "gene", "score", "strand"]]

# write sff, gff, and bed files
print('Writing chloe outputs')
sff.to_csv(args.output + '/result.sff', sep="\t", header=False, index=False)
gff.to_csv(args.output + '/result.gff', sep="\t", header=False, index=False)
bed.to_csv(args.output + '/result.bed', sep="\t", header=False, index=False)


# function to read fasta
def read_fasta(filename):
    name, seq = None, ""
    fasta = open(filename, 'r')
    for line in fasta:
        if line.startswith('>') and name is None:
            name = line.rstrip('\n').replace('>', '')
        else:
            if line.startswith('>') and name is not None:
                yield [name, seq]
                name = line.rstrip('\n').replace('>', '')
                seq = ''
            else:
                seq = seq + line.rstrip('\n')
    yield [name, seq]
    fasta.close()


print('Generating fasfa file')

# open bed to read
bed = open(args.output + '/result.bed', 'r')

# open fas to write
fas = open(args.output + '/result.fas', 'w')

# iterate through lines in the bed file
for line in bed:
    line = line.rstrip('\n').split('\t')
    sequence_id = line[0]
    start = int(float(line[1]))
    end = int(float(line[2]))
    gene = line[3]
    score = line[4]
    orientation = line[5]

    # write folded sequences to fasta
    for fas_name, fas_seq in read_fasta('NC_000932.1.fasta'):
        if sequence_id == fas_name.split(' ')[0]:
            fas.write('>' + sequence_id + '; ' + str(start) + "-" + str(end) + '; ' + orientation + '; ' + gene + '\n')
            tmp_seq = fas_seq[start:end]
            fold_len = 60
            for i in range(0, len(tmp_seq), fold_len):
                fas.write(tmp_seq[i:i + fold_len] + '\n')

# close files
bed.close()
fas.close()

