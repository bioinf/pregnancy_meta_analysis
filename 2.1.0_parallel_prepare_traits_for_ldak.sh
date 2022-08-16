#!/bin/bash

DATENOW=$( date )
echo "Started at: $DATENOW"

THREADS=$1
for ((i=0;i<=THREADS;i++)); do
  python ./2.0.0_prepare_traits_for_ldak.py ${i} ${THREADS} &
done
wait

DATENOW=$( date )
echo "Finished at: $DATENOW"
