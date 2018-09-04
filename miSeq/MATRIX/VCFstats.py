#! /bin/env python

# From VEP annotated VCF, it extracts the info of interest and outputs a .csv with corresponding columns

import sys
import re
import argparse
import pandas as pd
import numpy as np
parser=argparse.ArgumentParser(description='Arguments to run VCFstats.py')
parser.add_argument('--file', required=True)
parser.add_argument('--sample', required=True)
args=parser.parse_args()


# Open VEP annotated VCF
filename=args.file
infile=open(filename,'r')
basename=filename.split('.')[0]

# Prepare output
# Create the header of the output .csv ### Adding Sample name in the first column ('SAMPLE')
OUTPUT_HEADER = ['SAMPLE','CHROM','POS','REF','ALT','FILTER','DP','GT','AF','Allele','Consequence','IMPACT','SYMBOL','Gene','Feature_type','Feature','BIOTYPE','VARIANT_CLASS','Protein_position','Amino_acids','Codons','Existing_variation']#,'gnomAD_AF','gnomAD_AFR_AF','gnomAD_AMR_AF','gnomAD_ASJ_AF','gnomAD_EAS_AF','gnomAD_FIN_AF','gnomAD_NFE_AF','gnomAD_OTH_AF','gnomAD_SAS_AF','CLIN_SIG','SOMATIC','PHENO'] # check the INFO header of your VCF for more ouput options
# Dataframe to store and manage all the data
dtf = pd.DataFrame(columns = OUTPUT_HEADER)
SAMPLE_name=args.sample

for l in infile:

    if l.startswith('##INFO=<ID=CSQ,'): # get the info form the CSQ header
        CSQ_HEADER=re.compile('Format: (.*)">').search(l).group(1).split('|')
        continue
    if l.startswith('##'): continue     #other '##' lines unused
    if l.startswith('#CHROM'): #'#CHROM' is the header
        HEADER=l.strip('#').split()
        SAMPLES=HEADER[9:] # The sample name is on the 9th (and next) field of the VCF file header.
        continue

    s=l.strip().split('\t')  # s will contain all fields for this line (including CSQ)
    s=dict(zip(HEADER, s)) # Store all fields about a variant in a dictionary
    for k in SAMPLES: s[k]=s[k] # Add a new field with the sample name

    # Extract sample info of interest: DP, GT, AF (CHANGE IT IF MORE INFO IS NEEDED). IMPORTANT: Assuming only one sample (working with a single sample in the VCF
    ## CSQ
    SAMPLE = SAMPLES[0] # using the first sample. This script assumes only one sample (the tumor) in the input VEP annotated VCF file. Consider modifications if needed
    SAMPLE_INFO = dict(zip(s['FORMAT'].split(':'),s[SAMPLE].split(':'))) # Ex: {'GT': '0/1', 'AD': '33,2', 'AF': '0.081', ... , 'SA_MAP_AF': '0.00,0.061,0.057', 'SA_POST_PROB': '0.011,0.019,0.970'}
    csq_aux = [tuple(x.split('=')) for x in s['INFO'].split(';')]  # Ex: [('DP', '37'), ('ECNT', '1'), ('POP_AF', '5.000e-08'), ('P_CONTAM', '0.00'), ('P_GERMLINE', '-8.142e+00'), ('TLOD', '4.21'), ('CSQ', 'T|intron_v...]
    csq_aux = [y for y in csq_aux if len(y)>1]
    SAMPLE_INFO_CSQ = dict(csq_aux)['CSQ'].split(',')   # Ex: ['T|intron_variant|MODIFIER|STAG2|ENSG00000101972|Transcript|ENST00000218089|p...]

    # Parameters of interest to be included in the output file
    DP = dict(csq_aux)['DP'] # sequencing depth (sa DP in VCF info)
    GT = SAMPLE_INFO['GT']
    AF = SAMPLE_INFO['AF']

    if len(SAMPLE_INFO_CSQ) > 1:
        CSQ = dict(zip(CSQ_HEADER,SAMPLE_INFO_CSQ[0].split('|'))) # Only retrieves the first VEP annotation
    else:
        CSQ = dict(zip(CSQ_HEADER,SAMPLE_INFO_CSQ[0].split('|')))

    # Finally assemble all the info to be written for each variant:
    OUTPUT = [SAMPLE_name, s['CHROM'],s['POS'],s['REF'],s['ALT'],s['FILTER'],DP,GT,AF]
    CSQ= [CSQ['Allele'],CSQ['Consequence'],CSQ['IMPACT'],CSQ['SYMBOL'],CSQ['Gene'],CSQ['Feature_type'],CSQ['Feature'],CSQ['BIOTYPE'],CSQ['VARIANT_CLASS'],CSQ['Protein_position'],CSQ['Amino_acids'],CSQ['Codons'],CSQ['Existing_variation']] #,CSQ['gnomAD_AF'],CSQ['gnomAD_AFR_AF'],CSQ['gnomAD_AMR_AF'],CSQ['gnomAD_ASJ_AF'],CSQ['gnomAD_EAS_AF'],CSQ['gnomAD_FIN_AF'],CSQ['gnomAD_NFE_AF'],CSQ['gnomAD_OTH_AF'],CSQ['gnomAD_SAS_AF'],CSQ['CLIN_SIG'],CSQ['SOMATIC'],CSQ['PHENO']] # Create a list with the fiels from CSQ to be written in the output
    OUTPUT = OUTPUT + CSQ
    dtf.loc[len(dtf)] = OUTPUT

#####################

# Save to a new file (.csv) which will contain the FILTERED results
dtf.to_csv('-'.join(['/home/jlanillos/Disco4tb/Projects/miSeq/results/VCF/SAMPLES2/allVariants.csv']),sep='\t', mode = 'a')





#df = pd.read_csv(filename, sep='\t')
#names = list(df.SAMPLE.unique())

#only = pd.DataFrame(index = names)
#categ = df.FILTER.unique()

#df.groupby(['SAMPLE', 'FILTER']).size()
#df['FILTER'].value_counts()
#t = df.loc[df['FILTER'] == 'clustered_events']


#PASS = df.loc[df['FILTER'] == 'PASS']
#CLUS = df.loc[df['FILTER'] == 'clustered_events']
#TLOD = df.loc[df['FILTER'] == 't_lod']

#AF_pass = pd.to_numeric(PASS['AF'])
#AF_clus = pd.to_numeric(CLUS['AF'])
#AF_tlod = pd.to_numeric(TLOD['AF'])






#In [84]: AF_tlod.mean()
#Out[84]: 0.0788163620807666

#In [85]: AF_clus.mean()
#Out[85]: 0.21486030871145503

#In [86]: AF_pass.mean()
#Out[86]: 0.5130137673994066
