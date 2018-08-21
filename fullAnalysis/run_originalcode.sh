#!/bin/bash

if [ "$#" -lt 2 ]; then
  echo "Usage: ./run_originalcode.sh <datapath> <dataset>"
  echo "Ex: ./run_originalcode.sh  /my/path/to/data DS1"
  exit 1
fi

datapath="$1"
TOTEMdataset="$2"

OUTPUT=fulloriginal-times.csv
headers="Date,Test,Threads,Execution time"
nthreads=0

echo $headers >> $OUTPUT

# Diagonal is fixed for this study
diagonal=-d45b

# Clone repo
git clone https://github.com/JavierCvilla/analysis_elastic.6500GeV.beta90.10sigma.git
cd analysis_elastic.6500GeV.beta90.10sigma/
git checkout remove-path

# Adjust Path to data
find -type f -name "input*h" | xargs sed -i "s#{DATAPATH}#$datapath#"

# Path to generators.root
sed -i "s#{PWD}#$PWD#" common_algorithms.h

# Set storage dir
sed -i "s#{PWD}#$PWD#" $TOTEMdataset/block1/parameters.h

# Setup environment
source /cvmfs/sft.cern.ch/lcg/views/dev3/latest/x86_64-centos7-gcc62-opt/setup.sh

# Run distill
date=$(date +%d-%m-%Y\ %H:%M)
echo "Running distill.cc: ./run $diagonal distill.cc $TOTEMdataset/block1/"
./run distill.cc $diagonal -no-bckg $TOTEMdataset/block1/
realtime=`cat $TOTEMdataset/block1/distill_45b_56t.log_run | grep real | awk '{print $NF}' | tr "ms" ": "`

echo "$date,distill,$nthreads,$realtime" >> $OUTPUT

# Run distributions
echo "Running distributions.cc"
date=$(date +%d-%m-%Y\ %H:%M)
./run distributions.cc $diagonal -no-bckg dTOTEMdataset/block1/
realtime=`cat $TOTEMdataset/block1/distributions_45b_56t.log_run | grep real | awk '{print $NF}' | tr "ms" ": "`

echo "$date,distill,$nthreads,$realtime" >> $OUTPUT

