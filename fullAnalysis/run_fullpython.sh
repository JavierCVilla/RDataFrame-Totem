#!/bin/bash

SCRIPTPATH=$(dirname "$0")
INPUTPATH=`dirname $SCRIPTPATH`

OUTPUT="fullpython-times.csv"

headers="Date,Dataset,Test,Threads,Execution time"

echo $headers >> $OUTPUT

for nthreads in 0 1 2 4 8 16
do
  th_tag="no_MT"
  if [ $nthreads -gt 0 ]; then
    th_tag="threads_$nthreads"
  fi

  # Distill
  echo "Running: distill.py - $th_tag"
  date=$(date +%d-%m-%Y\ %H:%M)
  # Store real execution time and discard binary output
  logfile="distill_DS1_${th_tag}_d45b_56t_new.log"
  realtime=$((/usr/bin/time -f'%E' python distill.py d45b_56t DS1 ${nthreads//0/" "} >> $logfile ) 2>&1 )
  # Supress warning during snapshot
  realtime=`echo -e $realtime | tail -n 1 | awk '{print $NF}'`
  echo "$date,DS1,distill,$nthreads,$realtime" >> $OUTPUT

  # Distributions
  echo "Running: distributions.py"

  INPUTFILE="distill_DS1_${th_tag}_d45b_56t_new.root"
  realtime=$((/usr/bin/time -f'%E' python distributions.py $INPUTFILE $nthreads  > /dev/null )     2>&1 )
  echo "$date,DS1,distributions,$nthreads,$realtime" >> $OUTPUT
  
done

