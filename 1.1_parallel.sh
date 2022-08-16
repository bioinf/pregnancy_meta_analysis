#!/bin/bash

DATENOW=$( date )
echo "Started at: $DATENOW"


MAP_FILE=$1
TEMPLATE=$2

while IFS="," read -r code ukb finn n
do
  while [ $( ps -f -u $USER | grep 'pipe_go' | wc -l ) -ge 3 ]; do sleep 1; done
  DIR=$( pwd; )/${TEMPLATE}_$code
  echo $DIR
  mkdir -p ${TEMPLATE}_$code
  cd ${TEMPLATE}_$code
  mkdir -p data images scripts
  python3 ../1.0.0_pipe_go_R6.py \
              $ukb \
              $finn \
              ${DIR}/data \
              ${DIR}/images \
              $n  \
              ${DIR}/scripts \
              $code \
             'SAMPLESIZE' \
             'false' \
             'false' > log.txt &
  cd ..
  
done < <(cat $MAP_FILE)
wait

DATENOW=$( date )
echo "Finished at: $DATENOW"
