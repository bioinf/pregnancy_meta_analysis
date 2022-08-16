#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import os

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
    "z_score",
]


VARIANTS = "/media/array/pregnancy/data/variants.tsv"


# In[2]:


variants = pd.read_csv(VARIANTS, sep='\t')


# In[3]:


dirs = [f for f in os.listdir(".") if os.path.isdir(f) and "analysis_n" in f]


# In[4]:

for cur_dir in dirs:
    print(f"Processing {cur_dir}")
    MTAG_FILE = f"{cur_dir}/data/_mtag_meta.txt"
    EXT_MTAG_FILE = f"{cur_dir}/data/ext_mtag_meta.txt"

    met_data = pd.read_csv(MTAG_FILE, sep='\t')

    met_merged_data = met_data.merge(variants[['rsid', 'chr', 'pos', 'ref', 'alt']], left_on='SNP', right_on='rsid')

    met_map_cols = {'mtag_z':'z_score','mtag_pval':'pval', "mtag_beta":"beta", "mtag_se": "se", "meta_freq": "maf"}
    met_merged_data = met_merged_data.rename(columns=met_map_cols)
    met_merged_data.maf = met_merged_data.maf.apply(lambda x: min(x, 1-x))
    print('\nSaving file...')
    met_merged_data[needed_cols].to_csv(EXT_MTAG_FILE, index=False, sep='\t')
    print(f'\nSaved to {EXT_MTAG_FILE}.')

print('Pipeline successfully finished!\n-----------------------------\n')
