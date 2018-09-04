# for i in vep*; do bgzip $i; tabix -p vcf $i.gz;done

# In fact, merging must be done with VEP annotated VCF files
# bcftools merge -o MERGED_all.vcf -O v *    OR            bcftools merge -i DP:join -o MERGED_all.vcf -O v vep1.vcf.gz vep2.vcf.gz
# bgzip -c MERGED_all.vcf > MERGED_all.vcf.gz;tabix -p vcf MERGED_all.vcf.gz; bcftools norm -m -both -o MERGED_multiallelic_all.vcf.gz -Oz MERGED_all.vcf.gz

# Read the merged_all VCF (contains all variants form al lsamples merged using bcftools)
parser=argparse.ArgumentParser(description='Arguments to processVEP.py')
parser.add_argument('--file', required=True)
args=parser.parse_args()
filename=args.file
infile=open(filename,'r')
basename=filename.split('.')[0]

import sys
import re
import argparse
import pandas as pd
import numpy as np



for l in infile:

#metadata
    if l.startswith('##INFO=<ID=CSQ,'): # get the info form the CSQ header
        CSQ_HEADER=re.compile('Format: (.*)">').search(l).group(1).split('|')
        continue
    if l.startswith('##'): continue     #other '##' lines unused
    if l.startswith('#CHROM'): #'#CHROM' is the header
        HEADER=l.strip('#').split()
        SAMPLES=HEADER[9:] # The sample name is on the 9th (and next) field of the VCF file header.
            # Create the Output HEADER (Filter: it's the one established after merging VCF files with bcftools merge (default); DP_allsamples (sum of the DP of all samples containing that variant (bftools merge -i DP:sum)
        OUTPUT_HEADER = ['ID','CHROM','POS','REF','ALT','FILTER','DP_allsamples','Allele','Consequence','IMPACT','SYMBOL','Gene','Feature_type','Feature','BIOTYPE','EXON','cDNA_position','VARIANT_CLASS','HGNC_ID','Protein_position','Amino_acids','Codons','Existing_variation','gnomAD_AF','gnomAD_AFR_AF','gnomAD_AMR_AF','gnomAD_ASJ_AF','gnomAD_EAS_AF','gnomAD_FIN_AF','gnomAD_NFE_AF','gnomAD_OTH_AF','gnomAD_SAS_AF','CLIN_SIG','SOMATIC','PHENO'] # check the INFO header of your VCF for more ouput options
        SAMPLES_HEADER = []
        for sample in SAMPLES: SAMPLES_HEADER = SAMPLES_HEADER + ['_'.join([sample,'gt'])] + ['_'.join([sample,'ad'])] + ['_'.join([sample,'af'])]
        AF_samplenames = []
        for sample in SAMPLES: AF_samplenames = AF_samplenames + ['_'.join([sample,'af'])] # SAMPLE_NAMES_af to be used later during filtering
        OUTPUT_HEADER = OUTPUT_HEADER + SAMPLES_HEADER
        OUTPUT= []
        #df = pd.DataFrame(columns = OUTPUT_HEADER)
        continue

#variants

    # s will contain all fields for this line (including CSQ); Store all fields of a variant in a dictionary
    s=l.strip().split('\t')
    s=dict(zip(HEADER, s))

    # Add a new field containing each sample's name
    for k in SAMPLES: s[k]=s[k]

# Extract variant fields
    ## CSQ field
    csq_aux = [tuple(x.split('=')) for x in s['INFO'].split(';')]  # Ex: [('DP', '37'), ('ECNT', '1'), ('POP_AF', '5.000e-08'), ('P_CONTAM', '0.00'), ('P_GERMLINE', '-8.142e+00'), ('TLOD', '4.21'), ('CSQ', 'T|intron_v...]
    csq_aux = [y for y in csq_aux if len(y)>1]
    # SAMPLE_INFO_DP: the addition of all samples' sequencing depth
    SAMPLE_INFO_DP = int(dict(csq_aux)['DP'].split(',')[0])
    # SAMPLE_INFO_CSQ. useful to extract transcripts info. Ex: ['T|intron_variant|MODIFIER|STAG2|ENSG00000101972|Transcript|ENST00000218089|p...]
    SAMPLE_INFO_CSQ = dict(csq_aux)['CSQ'].split(',')
    # Parameters of interest to be included in the output file
    ID = '_'.join([s['CHROM'],s['POS'],s['REF'],s['ALT']]) # Variant ID containing the chrom, pos, ref and alt


# Extract information per sample
    SAMPLE_INFO_fields = []
    OUTPUT_general = [ID, s['CHROM'],s['POS'],s['REF'],s['ALT'],s['FILTER'],SAMPLE_INFO_DP]
    for SAMPLE in SAMPLES:
        SAMPLE_INFO = dict(zip(s['FORMAT'].split(':'),s[SAMPLE].split(':'))) # Ex: {'GT': '0/1', 'AD': '33,2', 'AF': '0.081', ... , 'SA_MAP_AF': '0.00,0.061,0.057', 'SA_POST_PROB': '0.011,0.019,0.970'}
        SAMPLE_INFO_fields = SAMPLE_INFO_fields + [SAMPLE_INFO['GT'], SAMPLE_INFO['AD'], SAMPLE_INFO['AF']] # SAMPLE_INFO['AD'][0], SAMPLE_INFO['AD'][1] are ref and alt allelic depths, respectively

# Preparing the output line and reading each annotated transcript for each variant. Sample info from CSQ added (see for loop below)
    for transcript in SAMPLE_INFO_CSQ: # There might be more than one transcript annotated, so they will be written in different lines
        CSQ = dict(zip(CSQ_HEADER,transcript.split('|')))
        CSQ= [CSQ['Allele'],CSQ['Consequence'],CSQ['IMPACT'],CSQ['SYMBOL'],CSQ['Gene'],CSQ['Feature_type'],CSQ['Feature'],CSQ['BIOTYPE'],CSQ['EXON'],CSQ['cDNA_position'],CSQ['VARIANT_CLASS'],CSQ['HGNC_ID'],CSQ['Protein_position'],CSQ['Amino_acids'],CSQ['Codons'],CSQ['Existing_variation'],CSQ['gnomAD_AF'],CSQ['gnomAD_AFR_AF'],CSQ['gnomAD_AMR_AF'],CSQ['gnomAD_ASJ_AF'],CSQ['gnomAD_EAS_AF'],CSQ['gnomAD_FIN_AF'],CSQ['gnomAD_NFE_AF'],CSQ['gnomAD_OTH_AF'],CSQ['gnomAD_SAS_AF'],CSQ['CLIN_SIG'],CSQ['SOMATIC'],CSQ['PHENO']] # Create a list with the fiels from CSQ to be written in the output

        OUTPUT.append(OUTPUT_general + CSQ + SAMPLE_INFO_fields)

# OUTPUT list stored into pandas df
df = pd.DataFrame(OUTPUT,columns = OUTPUT_HEADER)

######## FILTERING THE DATAFRAME

    # FILTER (PASS and Clustered_events pass the filter)
#filt = df.loc[(df['FILTER'] == 'PASS') | (df['FILTER'] == 'clustered_events')] # (is 'PASS' OR 'clustered_events')--> the problem of using this line is that the filter value assigned to each
# variant gives information about the first variant annotated..so, if there are two samples with the same variant, the filter filed is annotated with the info of the first sample with that variant found.


    # Filtering by Allele Frequency (AF) > 0.15. Preliminary steps: convert '.' to NaN --> convert AF columns into float --> Apply AF filtering
for sampleAF in AF_samplenames: filt[[sampleAF]] = filt[[sampleAF]].replace('.',np.nan)
for sampleAF in AF_samplenames: filt[sampleAF] = pd.to_numeric(filt[sampleAF]).values # Allele frequencies column (str) to float

filt = filt.loc[(filt[AF_samplenames] >= 0.15).any(axis=1)]

#for sampleAF in AF_samplenames: filt = filt.loc[(filt[sampleAF]  >= 0.15) | (np.isnan(filt[sampleAF])) ] # >> IMPORTANT: This line serves only if you are interested on filtering variants with AF> i nall the samples containing such variant
# Filter by consequence ('missense_variant','frameshift_variant','stop_gained','missense_variant&splice_region_variant','missense_variant&NMD_transcript_variant',
    # 'stop_gained&NMD_transcript_variant','inframe_deletion','stop_gained&frameshift_variant','stop_gained&frameshift_variant&NMD_transcript_variant',
    # 'missense_variant&splice_region_variant&NMD_transcript_variant','stop_lost','inframe_insertion','inframe_insertion&NMD_transcript_variant'
filt = filt[filt['Consequence'].str.contains('missense_variant|frameshift_variant|stop_gained|missense_variant&splice_region_variant|missense_variant&NMD_transcript_variant|stop_gained&NMD_transcript_variant|inframe_deletion|stop_gained&frameshift_variant|stop_gained&frameshift_variant&NMD_transcript_variant|missense_variant&splice_region_variant&NMD_transcript_variant|stop_lost|inframe_insertion|inframe_insertion&NMD_transcript_variant')==True]

#t_lod
filt_tlod = df.loc[(df['FILTER'] == 't_lod')] # (is 't_lod')
for sampleAF in AF_samplenames: filt_tlod[[sampleAF]] = filt_tlod[[sampleAF]].replace('.',np.nan)
#filt_tlod = filt_tlod.loc[(filt_tlod[AF_samplenames] >= 0.15).any(axis=1)]
#for sampleAF in AF_samplenames: filt_tlod[sampleAF] = pd.to_numeric(filt_tlod[sampleAF]).values # Allele frequencies column (str) to float
#for sampleAF in AF_samplenames: filt_tlod = filt_tlod.loc[filt_tlod[sampleAF]  >= 0.1 | (np.isnan(filt[sampleAF]))]
filt_tlod = filt_tlod[filt_tlod['Consequence'].str.contains('missense_variant|frameshift_variant|stop_gained|missense_variant&splice_region_variant|missense_variant&NMD_transcript_variant|stop_gained&NMD_transcript_variant|inframe_deletion|stop_gained&frameshift_variant|stop_gained&frameshift_variant&NMD_transcript_variant|missense_variant&splice_region_variant&NMD_transcript_variant|stop_lost|inframe_insertion|inframe_insertion&NMD_transcript_variant')==True]


##################### APPRIS #####################

# Open APPRIS database
apprisdb = pd.read_csv('/home/jlanillos/Disco4tb/refGen/APPRIS/appris_data.principal.txt', header=None, sep='\t')
# Find transcripts in our results that are stored in filt = filt.loc[(filt[AF_samplenames] >= 0.15).any(axis=1)]APPRIS_data.principal database
# Ángel's way
apprisdictPrincipal = dict(zip(apprisdb[2],apprisdb[4]))
apprisdictCCDS = dict(zip(apprisdb[2],apprisdb[3]))
# For loop to modify the transcript tags (e.g., PRINCIPAL:4 --> PRINCIPAL_4
for key in apprisdictPrincipal.keys():
    apprisdictPrincipal[key] = apprisdictPrincipal[key].replace(':','_')
# Add two columns containing transcript tags (like PRINCIPAL:4) and associated CCDS
filt['CCDS'] = ''
filt['APPRIS'] = ''
filt.CCDS = filt['Feature'].map(apprisdictCCDS)
filt.APPRIS = filt['Feature'].map(apprisdictPrincipal)

filt_tlod['CCDS'] = ''
filt_tlod['APPRIS'] = ''
filt_tlod.CCDS = filt['Feature'].map(apprisdictCCDS)
filt_tlod.APPRIS = filt['Feature'].map(apprisdictPrincipal)

df['CCDS'] = ''
df['APPRIS'] = ''
df.CCDS = filt['Feature'].map(apprisdictCCDS)
df.APPRIS = filt['Feature'].map(apprisdictPrincipal)

##################### SAVE RESULTS #####################

# outdir = '/home/jlanillos/Disco4tb//Projects/miSeq/src/MATRIX20180820/results/' # testing the code
outdir = '/home/jlanillos/Disco4tb/Projects/miSeq/results_SAMPLES3/CORRECTED/MATRIX_GT_AD_AF/'

# NON-FILTERED OUTPUT
df.to_csv(''.join([outdir,'UNFILTERED_MATRIX.csv']),sep='\t')

# FILTER TLOD: includes only t_lod (no AF filter applied)
filt_tlod.to_csv(''.join([outdir,'TLOD_MATRIX.csv']),sep='\t')

# FINAL FILTERED OUTPUT
filt.to_csv(''.join([outdir,'FINAL_MATRIX.csv']),sep='\t')





