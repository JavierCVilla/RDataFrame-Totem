#!/bin/sh

TEST="distributions-allhistos"
BINS=("$TEST-jit" "$TEST-jit-opt" "$TEST-nojit" "$TEST-nojit-opt")
SCRIPTPATH=$(dirname "$0")
INPUTPATH=`dirname $SCRIPTPATH`
INPUTFILE="$INPUTPATH/distill_DS1_d45b_56t.root"

OUTPUT="$TEST-times.csv"

date=$(date +%d-%m-%Y\ %H:%M)
headers="Date,Test,Threads,Execution time"

echo "Execution: $date;" >> $OUTPUT
echo $headers >> $OUTPUT

for binary in "${BINS[@]}"
do
   for nthreads in 0 1 2 4 8
   do
     echo "Running: $binary - $nthreads threads" 
     date=$(date +%d-%m-%Y\ %H:%M)
     # Store real execution time and discard binary output
     realtime=$((/usr/bin/time -f'%E' ./$binary $INPUTFILE $nthreads  > /dev/null ) 2>&1 )
     echo "$date,$binary,$nthreads,$realtime" >> $OUTPUT
   done
done

