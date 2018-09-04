#! /bin/env python

# From VEP annotated VCF, it extracts the info of interest and outputs a .csv with corresponding columns

import sys
import re
import argparse
import pandas as pd
import numpy as np
parser=argparse.ArgumentParser(description='Arguments to processVEP.py')
parser.add_argument('--file', required=True)
args=parser.parse_args()


filename=args.file
infile=open(filename,'r')
basename=filename.split('.')[0]

# Create the header of the output .csv
OUTPUT_HEADER = ['ID', 'CHROM','POS','REF','ALT','FILTER','DP','GT','AF','Allele','Consequence','IMPACT','SYMBOL','Gene','Feature_type','Feature','BIOTYPE','VARIANT_CLASS','Protein_position','Amino_acids','Codons','Existing_variation','gnomAD_AF','gnomAD_AFR_AF','gnomAD_AMR_AF','gnomAD_ASJ_AF','gnomAD_EAS_AF','gnomAD_FIN_AF','gnomAD_NFE_AF','gnomAD_OTH_AF','gnomAD_SAS_AF','CLIN_SIG','SOMATIC','PHENO'] # check the INFO header of your VCF for more ouput options
OUTPUT = []
# Dataframe to store and manage all the data
dtf = pd.DataFrame(columns = OUTPUT_HEADER)
infile=open(filename,'r')
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
    ID = '_'.join([s['CHROM'],s['POS'],s['REF'],s['ALT']])
    DP = dict(csq_aux)['DP'] # sequencing depth (sa DP in VCF info)
    GT = SAMPLE_INFO['GT']
    AF = SAMPLE_INFO['AF']

    if len(SAMPLE_INFO_CSQ) > 1:
        CSQ = dict(zip(CSQ_HEADER,SAMPLE_INFO_CSQ[0].split('|'))) # Only retrieves the first VEP annotation
    else:
        CSQ = dict(zip(CSQ_HEADER,SAMPLE_INFO_CSQ[0].split('|')))

    # Finally assemble all the info to be written for each variant:
    OUTPUT1 = [ID, s['CHROM'],s['POS'],s['REF'],s['ALT'],s['FILTER'],DP,GT,AF]
    CSQ= [CSQ['Allele'],CSQ['Consequence'],CSQ['IMPACT'],CSQ['SYMBOL'],CSQ['Gene'],CSQ['Feature_type'],CSQ['Feature'],CSQ['BIOTYPE'],CSQ['VARIANT_CLASS'],CSQ['Protein_position'],CSQ['Amino_acids'],CSQ['Codons'],CSQ['Existing_variation'],CSQ['gnomAD_AF'],CSQ['gnomAD_AFR_AF'],CSQ['gnomAD_AMR_AF'],CSQ['gnomAD_ASJ_AF'],CSQ['gnomAD_EAS_AF'],CSQ['gnomAD_FIN_AF'],CSQ['gnomAD_NFE_AF'],CSQ['gnomAD_OTH_AF'],CSQ['gnomAD_SAS_AF'],CSQ['CLIN_SIG'],CSQ['SOMATIC'],CSQ['PHENO']] # Create a list with the fiels from CSQ to be written in the output
    OUTPUT1 = OUTPUT1 + CSQ
    OUTPUT.append(OUTPUT1)



    #dtf.loc[len(dtf)] = OUTPUT
dtf = pd.DataFrame(OUTPUT,columns = OUTPUT_HEADER)



dtf.to_csv('-'.join([basename,'annotations.csv']),sep='\t')


### Create another annotation file with some filters (AF, mutation type, variant caller filter...)

##### FILTERING #####
dtf['AF'] = pd.to_numeric(dtf.AF).values # Allele frequencies column (str) to float


filt = dtf.loc[dtf['FILTER'] == 'PASS']
filt['AF'] = pd.to_numeric(filt.AF).values # Allele frequencies column (str) to float


filt1 = filt.loc[filt['AF'] >= 0.1]
#filt = filt.loc[filt['Consequence'] != 'intron_variant']

#filt = filt.loc[(filt['Consequence'] != 'intron_variant') & (filt['Consequence'] != '5_prime_UTR_variant') ]
filt2 = filt1[filt1['Consequence'].str.contains('missense_variant|frameshift_variant|stop_gained|missense_variant&splice_region_variant|missense_variant&NMD_transcript_variant|stop_gained&NMD_transcript_variant|inframe_deletion|stop_gained&frameshift_variant|stop_gained&frameshift_variant&NMD_transcript_variant|missense_variant&splice_region_variant&NMD_transcript_variant|stop_lost|inframe_insertion|inframe_insertion&NMD_transcript_variant')==True]

filt3  = dtf.loc[dtf['AF'] >= 0.1]

filt4 = filt3[filt3['Consequence'].str.contains('missense_variant|frameshift_variant|stop_gained|missense_variant&splice_region_variant|missense_variant&NMD_transcript_variant|stop_gained&NMD_transcript_variant|inframe_deletion|stop_gained&frameshift_variant|stop_gained&frameshift_variant&NMD_transcript_variant|missense_variant&splice_region_variant&NMD_transcript_variant|stop_lost|inframe_insertion|inframe_insertion&NMD_transcript_variant')==True]
#####################

# Save filtered results to a new file (.csv) which will contain the FILTERED results

filt.to_csv('-'.join([basename,'PASS_FILTERED_annotations.csv']),sep='\t')
filt1.to_csv('-'.join([basename,'AF_PASS_FILTERED_annotations.csv']),sep='\t')
filt2.to_csv('-'.join([basename,'ALL_FILTERED_annotations.csv']),sep='\t')
filt3.to_csv('-'.join([basename,'AF_FILTERED_annotations.csv']),sep='\t')
filt4.to_csv('-'.join([basename,'AF_CONSEQUENCE_FILTERED_annotations.csv']),sep='\t')



