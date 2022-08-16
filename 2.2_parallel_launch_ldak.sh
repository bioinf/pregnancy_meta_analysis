#!/bin/bash

DATENOW=$( date )
echo "Started at: $DATENOW"
i=0
# echo "cat"
for d in analysis_n_*/ ; do
    while [ $( ps -f -u $USER | grep 'ldak5.2' | wc -l ) -ge 27 ]; do sleep 1;  done
    echo -e "\n\n=============\nSTEP: ${i}        \n=============\n"
    i=$((i+1))
    echo ${d}
    cd ${d}
    ./ldak_script.sh &
    cd /media/MIRROR/ukb_finngen/meta_analysis/
done
wait

DATENOW=$( date )
echo "Finished at: $DATENOW"

