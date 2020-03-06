#!/bin/bash
if [ "$#" -ne 1 ]
then
  echo "Usage: $0 matrixName"
  exit 1
fi
matrix=$1

cd buildGnu
make clean
make -j

../runTest.sh parallelSpmvSm    $matrix gnu

cd ../

source setIcc intel64
source setImpi

cd buildIntel
make clean
make -j

../runTest.sh parallelSpmvSm    $matrix intel


cd ../

source setPgi 19.x
source setPgiMpi 19.x

cd buildPgi
make clean
make -j

../runTest.sh parallelSpmvSm    $matrix pgi

cd ../

