#!/usr/bin/env python
# coding: utf-8

# In[151]:


import pandas as pd
import numpy as np
import os
from tqdm import tqdm
from sys import argv


# In[163]:


VARIANTS = "/home/achangalidi/ukb_finn_extra_data/variants.tsv"
DATA_PATH_N = "/media/MIRROR/ukb_finngen/meta_analysis/file_mapping.csv"
DATA_PATH_N_PREG = "/media/MIRROR/ukb_finngen/meta_analysis/file_mapping_preg.csv"
TO_RENAME = {'Predictor':'Predictor', 'Allele1':'A1','Allele2':'A2','n':'n', 'Zscore':'Z'}

PREG_DIR = "/media/MIRROR/ukb_finngen/meta_analysis/PREG_DATA/"
CUR_DIR = "/media/MIRROR/ukb_finngen/meta_analysis/"

SCRIPT = '''/home/achangalidi/tools/ldak5.2.linux --sum-cors gencor \\
--summary {} \\
--summary2 {} \\
--tagfile /home/achangalidi/tools/ldak.thin.hapmap.gbr.tagging \\
--allow-ambiguous YES \\
--check-sums NO \\
> {}'''

CUR_PART = int(argv[1])
ALL_PARTS = int(argv[2])


# In[99]:


# variants = pd.read_csv(VARIANTS, sep='\t')
# path_n = pd.read_csv(DATA_PATH_N, header=None)
# path_n_preg = pd.read_csv(DATA_PATH_N_PREG, header=None)
# variants.head()


# In[142]:


## for every
def return_ldak_df(FILE):
    pair = FILE.split('/')[-1].replace('1.TBL.gz', '').replace('1.TBL', '')
    if len(path_n[path_n.iloc[:,0]==pair][3]) > 0:
        N = int(path_n[path_n.iloc[:,0]==pair][3])
    else:
        N = int(path_n_preg[path_n_preg.iloc[:,0]==pair][3])
        

    data = pd.read_csv(FILE, sep='\t')
    data_m = pd.merge(data, variants, left_on='MarkerName', right_on='rsid')
    data_m['Predictor'] = data_m.chr.astype(str)+':'+data_m.pos.astype(str)
    data_m['n'] = data_m['n_called']+N
    f_data = data_m.rename(columns = TO_RENAME)[list(TO_RENAME.values())]
    f_data = f_data[(f_data.A1.apply(lambda x: len(str(x))<2)) & (f_data.A2.apply(lambda x: len(str(x))<2))]
    f_data.A1 = f_data.A1.apply(lambda x: x.upper())
    f_data.A2 = f_data.A2.apply(lambda x: x.upper())
    return f_data


# In[143]:


FILES_PREG = [f for f in os.listdir(path=PREG_DIR) if f[-4:]=='.TBL']
preg_names = []
for f in FILES_PREG:
#     preg_data = return_ldak_df(os.path.join(PREG_DIR, f))
    name = f.replace('.TBL', '.tsv')
    out_file = os.path.join(PREG_DIR, name)
    print(out_file )
#     preg_data.to_csv(out_file , index=False, sep='\t')
    preg_names.append(out_file )


# In[182]:


ALL_DIRS = [f for f in os.listdir(path=CUR_DIR) if 'analysis_n_' in f]
batch_size = len(ALL_DIRS) // ALL_PARTS
ALL_DIRS = ALL_DIRS[batch_size*CUR_PART:batch_size*(CUR_PART+1)]

for cur_dir in tqdm(ALL_DIRS):
    cur_file = [f for f in os.listdir(path=os.path.join(cur_dir, "data")) if f[-7:]=='.TBL.gz' or f[-4:]=='.TBL'][0]
#     data2 = return_ldak_df(os.path.join(CUR_DIR, cur_dir, "data", cur_file))
##     data2.head()
    out_file = os.path.join(CUR_DIR, cur_dir, "data", cur_file.replace('.TBL', '.tsv').replace('.gz', ''))
    print(out_file)
#     data2.to_csv(out_file , index=False, sep='\t')


    script = f"WD=$(pwd)\n"
    script += f"cd {os.path.join(CUR_DIR, cur_dir, 'data')}\n\n\n\n"
    for preg_name in preg_names:
        script += SCRIPT.format(out_file, preg_name, f"{preg_name}_log.txt\n\n")
        script += f"rename 's/gencor/{preg_name.split('/')[-1].replace('.tsv', '').replace('.gz', '')}/g' * \n\n\n\n"
    script += "cd ${WD}\n"
    script_file = os.path.join(CUR_DIR, cur_dir, "ldak_script.sh")
    with open(script_file, 'w') as f:
        print(script, file = f)
    os.system(f"chmod 744 {script_file}")


