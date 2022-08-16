#!/bin/bash

DATENOW=$( date )
echo "Started at: $DATENOW"

DRAW_R_SCRIPT='/media/array/pregnancy/pregnancy_meta_analysis/1.0.1_DRAW.R'

CUR_DIR=$(pwd)
for d in analysis_n_*/ ; do
    while [ $( ps -f -u $USER | grep 'DRAW.R' | wc -l ) -ge 3 ]; do sleep 1;  done
    # d="analysis_O15_PREG_ECTOP/"
    cd ${d}
    # get current trait code
    CODE="${d/analysis_n_/""}"  
    CODE="${CODE/\//""}"
    echo ${CODE}
    # get finngen file
    pattern="${CUR_DIR}/${d}data/maf*hg19lifted*"
    finn_file=( $pattern )
    echo "${finn_file[0]}"
    
    
    MH_PLOT="Rectangular-Manhattan.pval_FinnG_${CODE}.pdf"
    QQ_PLOT="QQplot.pval_FinnG_${CODE}.pdf"
    
    # launch
    ${DRAW_R_SCRIPT} FinnG_${CODE} ${finn_file[0]} && mv ${MH_PLOT} images/ && mv ${QQ_PLOT} images/ && cd ${CUR_DIR} > log_fin.txt &
    cd ..
done
wait

DATENOW=$( date )
echo "Finished at: $DATENOW"


