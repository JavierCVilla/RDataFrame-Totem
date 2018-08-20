#!/bin/sh

TESTS=("all-histos" "one-histo" "no-histos")

echo "### Starting Test suite ###"
for t in "${TESTS[@]}"
do
  echo "- Running test: $t"
  cd $t
  export PATH=$PWD:$PATH
  runscript=run_${t/-/}.sh
  $runscript
  cd -
done
