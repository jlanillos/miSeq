#! /bin/env python

# From VEP annotated VCF, it extracts the info of interest and outputs a .csv with corresponding columns

import sys
import re
import argparse
import pandas as pd
import numpy as np
parser=argparse.ArgumentParser(description='Arguments to processVEP.py')
parser.add_argument('--file', required=True)
parser.add_argument('--out', required=True)
args=parser.parse_args()

out_dir = args.out #'/home/jlanillos/PycharmProjects/mTORTSC12/VCF/'
filename=args.file
dtf=pd.read_csv('/'.join([filename, 'vep_mutect-FILTERED_annotations.csv']),sep = '\t') # Input cvs

##### FILTERING #####

filt = dtf.loc[(dtf['SYMBOL'] == 'MTOR') | (dtf['SYMBOL'] == 'TSC1') | (dtf['SYMBOL'] == 'TSC2') ]
filt['SAMPLE'] = filename
cols = list(filt.columns.values)
cols = cols[-1:] + cols[:-1]
filt = filt[cols]

#####################

# Open a new file (.csv) which will append the results form all samples
filt.to_csv('/'.join([out_dir,'ANNOTATED_mTOR_TSC1_TSC2.csv']),sep='\t', mode = 'a')
