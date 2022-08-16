#!/usr/bin/env python3

import pandas as pd
from collections import Counter
import numpy as np
import os
from sys import argv


# UKB file with some additional info (rsids)
VARIANTS = "/media/array/pregnancy/data/variants.tsv"
# PATH to METAL
METAL_PATH = '/home/achangalidi/projects/pregnancy/generic-metal/executables/metal'
DRAW_R_SCRIPT = '/media/array/pregnancy/pregnancy_meta_analysis/1.0.1_DRAW.R'


# =========================
# "ARGS"
# =========================
# files with summstats

## It is critical that the ukb file go first, and then the finngen, 
## but this does not affect what happens in the METAL,
## because we preprocess first and then use the preprocessed files in METAL

NAME1 = argv[1]
NAME2 = argv[2]
# directories for output
OUT_DIR = argv[3]
IMAGES_DIR = argv[4]

# there is no column with amount of samples in Finngen (there is in meta data), 
# that is why we input it
DATA2_N_SAMPLES = int(argv[5])

MAF_CUT = 0.05

# file with METAL script for analysis (it will be self-generated)
METAL_SCRIPT_OUT = f'{argv[6]}/metal_launch_flc_fmaf1.txt'
# prefix of METAL output
METAL_OUT_BEGINNING = argv[7]

# how to run metal
HOW = argv[8]
if HOW in ['ss', 'SAMPLESIZE', 'samplesize', 'SS']:
    HOW = 'SAMPLESIZE'
elif HOW in ['se', 'STDERR', 'stderr', 'SE']:
    HOW = 'STDERR'
else:
    raise ValueError
    
    

# scale?? UKB betas and SEs
try:
    SCALE = argv[9]
    if SCALE in ['true', 'True', '1', 1]:
        print('scale')
        SCALE = True
    else:
        print('non-scale')
        SCALE = False
except IndexError:
    SCALE = False    
    

# delete evth unneded and gzip??
try:
    ECONOM = argv[10]
    if ECONOM in ['true', 'True', '1', 1]:
        print('eco')
        ECONOM = True
    else:
        print('non-eco')
        ECONOM = False
except IndexError:
    ECONOM = False




# =======================================================
# ADDITIONAL FILES and SETTINGS
# =======================================================
# Prepared summstats files for metal
FILE1_FOR_META = f'{OUT_DIR}/filtered_{NAME1.split("/")[-1]}_.tsv'
FILE2_FOR_META = f'{OUT_DIR}/filtered_{NAME2.split("/")[-1]}_.tsv'

FILE2_MAF = f'{OUT_DIR}/maf_fg_{NAME2.split("/")[-1]}_.tsv'

# METAL output 
METAL_OUT = f"{OUT_DIR}/{{}}{METAL_OUT_BEGINNING}1.TBL"
# for postprocessing METAL OUTPUT
BEFORE_PREFIX = ""
AFTER_PREFIX = "extended_"
# METAL output before postprocessing
METAL_FILE = METAL_OUT.format(BEFORE_PREFIX)
# METAL output after postprocessing
EXT_METAL_FILE = METAL_OUT.format(AFTER_PREFIX)




# =======================================================
# SECONDARY FUNCTIONS
# =======================================================

def delete_multialleles(data):
    def check(x):
        return x in multialleles

    count = pd.DataFrame.from_dict(Counter(data.rsid), orient="index").reset_index()
    multialleles = set(count[count[0] > 1]["index"])
    return data[~data.rsid.apply(check)]


def intersect_data(data1, data2):
    def needed_rsids(x):
        return x in common_rsids

    common_rsids = set(data1.rsid) & set(data2.rsid)
    print(f"rsids before intersection:\t{len(set(data1.rsid))}, {len(set(data2.rsid))}")
    print(f"rsids after intersection:\t{len(common_rsids)}")

    data1 = data1[data1.rsid.apply(needed_rsids)]
    data2 = data2[data2.rsid.apply(needed_rsids)]

    d1, d2 = set(data1.rsid), set(data2.rsid)
    print("Quick check: ", len(d1), len(d2), len(d1 & d2))
    print("Shapes: ", data1.shape, data2.shape)
    return data1, data2

def get_cols_to_save(data):
    try:
        data[additional_cols]
        return needed_cols + additional_cols
    except KeyError:
        return needed_cols
    
    
# For sorting purposes
chr_order = {'1': 0,
 '2': 1,
 '3': 2,
 '4': 3,
 '5': 4,
 '6': 5,
 '7': 6,
 '8': 7,
 '9': 8,
 '10': 9,
 '11': 10,
 '12': 11,
 '13': 12,
 '14': 13,
 '15': 14,
 '16': 15,
 '17': 16,
 '18': 17,
 '19': 18,
 '20': 19,
 '21': 20,
 '22': 21,
 'X': 22,
 'Y': 23}
    
   
    
# ============================
# Step0 unzipping
# ============================
unzip1_flag = False
print('Unzipping file1:')
print(NAME1, NAME2, sep='\n', end = '\n\n\n\n\n\n\n\n\n\n\n')
if NAME1[-4:] == '.bgz':
    NAME1_TO = os.path.join(OUT_DIR, NAME1.split('/')[-1][:-4])
    command = f'gunzip -c {NAME1} > {NAME1_TO}'
    print(command)
    os.system(command)
    NAME1 = NAME1_TO
    unzip1_flag = True
print('Done unzipping file_1!')

print('Unzipping file1:')
unzip2_flag = False
if NAME2[-3:] == '.gz':
    NAME2_TO = os.path.join(OUT_DIR, NAME2.split('/')[-1][:-3])
    command = f'gunzip -c {NAME2} > {NAME2_TO}'
    print(command)
    os.system(command)
    NAME2 = NAME2_TO
    unzip2_flag = True
print('Done unzipping file_2!')

  
    
    

# =======================================================
# STEP 1: Preprocessing
# Prepare ukb and finngen summstats for METAL format.
# =======================================================

needed_cols = [
    "rsid",
    "beta",
    "se",
    "pval",
    "maf",
    "chr",
    "pos",
    "ref",
    "alt",
    "n_sample",
    "z_score",
]

print('Reading datasets...')
data1 = pd.read_csv(NAME1, sep="\t")
data2 = pd.read_csv(NAME2, sep="\t")

print('\nProcessing dataset 1...')
data1["chr"] = data1.variant.apply(lambda string: string.split(":")[0])
data1["pos"] = data1.variant.apply(lambda string: str(int(string.split(":")[1])))
data1["ref"] = data1.variant.apply(lambda string: string.split(":")[2])
data1["alt"] = data1.variant.apply(lambda string: string.split(":")[3])

# here we get rsids from positions
variants = pd.read_csv(VARIANTS, sep="\t")
data1 = data1.merge(variants[["variant", "rsid"]], on="variant")
# ename columns to format
data1 = data1.rename(
    columns={
        "rsid": "rsid",
        "beta": "beta",
        "se": "se",
        "pval": "pval",
        "minor_AF": "maf",
        "#chrom": "chr",
        "pos": "pos",
        "ref": "ref",
        "alt": "alt",
        "n_complete_samples": "n_sample",
        "tstat": "z_score",
    }
)
# remove multi-alleles
data1 = delete_multialleles(data1)


print('\nProcessing dataset 2...')
data2["#chrom"] = data2["new_chr"].astype(str).apply(lambda x: x.replace('chr', ''))
data2["pos"] = data2["new_coord"].astype('float64').astype('Int64').astype(str)

# use hardcoded number of samples for finngen
data2["n_sample"] = DATA2_N_SAMPLES
# rename columns to format
data2["z_score"] = data2["beta"] / data2["sebeta"]
data2["af_alt"] = data2["af_alt"].apply(lambda x: min(x, 1-x))
data2 = data2.rename(
    columns={
        "rsids": "rsid",
        "beta": "beta",
        "sebeta": "se",
        "pval": "pval",
        "af_alt": "maf",
        "#chrom": "chr",
        "pos": "pos",
        "ref": "ref",
        "alt": "alt",
        "n_sample": "n_sample",
        "z_score": "z_score",
    }
)
# remove multi-alleles
data2 = delete_multialleles(data2)

print('\nRemoving nans, intersecting rsids...')
data1 = data1.dropna()
data2 = data2.dropna()


# sort by pos
data1['chr_order'] = data1["chr"].apply(lambda x: chr_order[x])
data1 = data1.sort_values(by=['chr_order', 'pos'])
data2 = data2[data2["chr"].isin(list(chr_order.keys()))]
data2['chr_order'] = data2["chr"].apply(lambda x: chr_order[x])
data2 = data2.sort_values(by=['chr_order', 'pos'])

# look for rs*,rs* in fg
data2['variant'] = data2.chr+':'+data2.pos+':'+data2.ref+':'+data2.alt
flags = data2.rsid.str.contains(',', case=False)
spec = data2[flags]
spec = spec.merge(variants[["variant", "rsid"]], on="variant", how = 'left')
spec = spec.dropna()
spec = spec.loc[spec.rsid_y.str.contains('rs', case=False)][['variant','rsid_y']]
data2 = data2.merge(spec, on='variant', how='left')
flag2 = ~data2.rsid_y.isna()
data2.loc[flag2, 'rsid'] = data2.loc[flag2, 'rsid_y']

# save the Finngen file only filtered by maf
data2[data2.maf >= MAF_CUT].loc[:,['chr', 'pos', 'ref', 'alt', 'rsid', 'pval', 'beta', 'se', 'maf']].to_csv(FILE2_MAF, index=False, sep='\t')
# ----

data1, data2 = intersect_data(data1, data2)

print('\nFiltering low confedent variants...')
data1_flc = data1[~data1.low_confidence_variant]
data1_flc, data2_flc = intersect_data(data1_flc, data2)

print(f'\nFiltering snps with maf < {MAF_CUT}')
data1_flc_fmaf1, data2_flc_fmaf1 = intersect_data(data1_flc[data1_flc.maf >= MAF_CUT], 
                                                  data2_flc[data2_flc.maf >= MAF_CUT])



if SCALE:
    print(f'\nTransforming beta and se of UKB data...')
    multiplier = data2_flc_fmaf1.se.std() / data1_flc_fmaf1.se.std()
    data1_flc_fmaf1.beta *= multiplier
    data1_flc_fmaf1.se *= multiplier

print('\nSaving data...')


data1_flc_fmaf1[needed_cols].to_csv(FILE1_FOR_META, index=False, sep='\t')
print(f'{FILE1_FOR_META} saved.')

data2_flc_fmaf1[needed_cols].to_csv(FILE2_FOR_META, index=False, sep='\t')
print(f'{FILE2_FOR_META} saved.')

print('\nDone!')




# =======================================================
# STEP 2: Launching METAL
# =======================================================


print('\nRunning METAL...')
# prepare a script that METAL will use]
metal_commands = f'''
SCHEME {HOW}
GENOMICCONTROL ON
AVERAGEFREQ ON
MINMAXFREQ ON

MARKER   rsid
WEIGHT   n_sample
ALLELE   ref alt
FREQ     maf
EFFECT   beta
STDERR   se
PVAL     pval

OUTFILE {METAL_OUT_BEGINNING} .TBL


PROCESS {FILE1_FOR_META}
PROCESS {FILE2_FOR_META}
# Execute meta-analysis
ANALYZE
'''


# цйму МЕТАД's script
with open(METAL_SCRIPT_OUT, 'w') as out_metal:
    out_metal.write(metal_commands)
    
# make it executable  
os.system(f'chmod 755 {METAL_SCRIPT_OUT}')
# execute it
os.system(f'{METAL_PATH} {METAL_SCRIPT_OUT}')
# METAL saves the file in the same place where it started, so let's move the resulting files to the data directory
os.system(f'mv  {METAL_OUT_BEGINNING}1.TBL {OUT_DIR}')
os.system(f'mv  {METAL_OUT_BEGINNING}1.TBL.info {OUT_DIR}')

print('\nMETAL done!')




# =======================================================
# STEP 3: Postprocessing METAL output
# =======================================================


needed_cols = [
    "rsid",
    "chr",
    "pos",
    "ref",
    "alt",
#     "z_score",
    "pval",
    #     "beta",
    #     "se",
    "maf",
    #     "n_sample",
]
additional_cols = [
    "beta",
    "se",
]


print('\nPostprocessing METAL out...')

met_data = pd.read_csv(METAL_FILE, sep='\t')
variants = pd.read_csv(VARIANTS, sep='\t')

# little game with columns to make it more convenient to use them later
met_merged_data = met_data.merge(variants[['rsid', 'chr', 'pos', 'ref', 'alt']], left_on='MarkerName', right_on='rsid')
met_map_cols = {'Zscore':'z_score','P-value':'pval', 'Freq1':'maf'}
met_merged_data = met_merged_data.rename(columns=met_map_cols)

met_merged_data.maf = np.min([met_merged_data.maf, 1-met_merged_data.maf], axis = 0)

print('\nSaving file...')
# take only the necessary columns and save
met_all_cols = get_cols_to_save(met_merged_data)
met_merged_data[met_all_cols].to_csv(EXT_METAL_FILE, index=False, sep='\t')

print(f'\nSaved to {EXT_METAL_FILE}.')



# =======================================================
# STEP 4: Drawing mh and qq plots.
# =======================================================



# just run the script on R (DRAW.R, lies nearby)
# for QQ and MH drawing
print(f'Drawing plots...')

os.system(f'{DRAW_R_SCRIPT} {METAL_OUT_BEGINNING} {EXT_METAL_FILE}')

MH_PLOT = f'Rectangular-Manhattan.pval_{METAL_OUT_BEGINNING}.pdf'
QQ_PLOT = f'QQplot.pval_{METAL_OUT_BEGINNING}.pdf'

os.system(f'mv {MH_PLOT} {IMAGES_DIR}')
os.system(f'mv {QQ_PLOT} {IMAGES_DIR}')




# =======================================================
# STEP 5: Deleting  shitty output
# =======================================================

    
if ECONOM:
    if unzip1_flag:
        os.system(f'rm {NAME1}')
    if unzip2_flag:
        os.system(f'rm {NAME2}')
    os.system(f'rm {FILE1_FOR_META}')
    os.system(f'rm {FILE2_FOR_META}')
    os.system(f'gzip {METAL_FILE}')
    print('gzipped')



print('Pipeline successfully finished!')