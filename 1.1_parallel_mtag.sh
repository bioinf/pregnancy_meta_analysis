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

  cd ${TEMPLATE}_$code

  python3 ../1.0.0_pipe_go_R6_mtag.py \
              $ukb \
              $finn \
              ${DIR}/data > log_mtag.txt &
  cd ..
  
done < <(cat $MAP_FILE)
wait

DATENOW=$( date )
echo "Finished at: $DATENOW"
