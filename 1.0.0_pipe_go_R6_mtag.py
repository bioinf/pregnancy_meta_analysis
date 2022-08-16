#!/usr/bin/env python3

from collections import Counter
import os
from sys import argv


# =========================
# "ARGS"
# =========================
NAME1 = "/media/array/pregnancy/data/HYPTENSPREG.gwas.imputed_v3.female.tsv"
NAME1 = argv[1]
NAME2 = "/media/array/pregnancy/data/finngen_R5_O15_HYPTENSPREG"
NAME2 = argv[2]
# directories for output
OUT_DIR = "/media/array/pregnancy/data"
OUT_DIR = argv[3]


# path to metal
MTAG_PATH = "/home/achangalidi/projects/pregnancy/mtag/mtag.py"


# =======================================================
# ADDITIONAL FILES and SETTINGS
# =======================================================
# Prepared summstats files for metal
FILE1_FOR_META = f'{OUT_DIR}/filtered_{NAME1.split("/")[-1]}_.tsv'
FILE2_FOR_META = f'{OUT_DIR}/filtered_{NAME2.split("/")[-1]}_.tsv'
 


# =======================================================
# STEP 2: Launching MTAG
# =======================================================

# запускаем метал
os.system(f'{MTAG_PATH} --sumstats  {FILE1_FOR_META},{FILE2_FOR_META} \
            --out {OUT_DIR}/ \
          --snp_name rsid \
          --beta_name beta \
          --se_name se \
          --p_name pval \
          --eaf_name maf \
          --chr_name chr \
          --bpos_name pos \
          --a1_name ref \
          --a2_name alt \
          --n_name n_sample \
          --z_name z_score \
          --force \
          --perfect_gencov \
          --equal_h2 \
          --no_overlap')

print('\nMETAL done!')