import argparse
import requests
import time
import pandas as pd

# argparse

argparse_description = '''
Python script to annotate chloroplast assemblies using Chloe https://chloe.plastid.org/
'''

argparse_usage = '''
python chloe_api.py --inpuut <input_fasta> --output <output_directory> 
'''

parser = argparse.ArgumentParser(description=argparse_description, usage=argparse_usage)
parser.add_argument("--input",  help="Input fasta to annotate",            required=True)
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
        names=["sequence_id", "source", "feature_type", "feature_start", "feature_end", "score", "strand", "phase", "attributes"])

# format bed file
bed = gff[ gff["feature_type"] == "gene" ]
bed["gene"] = bed["attributes"].str.split(";").str[1].str.replace("Name=","")
bed = bed[["sequence_id", "feature_start", "feature_end", "gene", "score", "strand"]]

# write files
print('Writing outputs')
sff.to_csv(args.output + '/result_sff.txt', sep="\t", header=False, index=False)
gff.to_csv(args.output + '/result_gff.txt', sep="\t", header=False, index=False)
bed.to_csv(args.output + '/result_bed.txt', sep="\t", header=False, index=False)

