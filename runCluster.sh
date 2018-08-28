#!/bin/bash
#   Crearte soft links with the programs into the same directory 
#   Do this for each compiler
#



if [ "$#" -ne 1 ]
then
  echo "Usage: $0 matrixName"
  exit 1
fi

cd buildGnu

../runClusterTest.sh parallelSpmv     $1 gnu
../runClusterTest.sh parallelSpmvSmV1 $1 gnu
../runClusterTest.sh parallelSpmvSmV2 $1 gnu

cd ../

source setIcc intel64
source setImpi

cd buildIntel

../runClusterTest.sh parallelSpmv     $1 intel
../runClusterTest.sh parallelSpmvSmV1 $1 intel
../runClusterTest.sh parallelSpmvSmV2 $1 intel


cd ../

source setPgi 18.x
source setPgiMpi 18.x

cd buildPgi

../runClusterTest.sh parallelSpmv     $1 pgi
../runClusterTest.sh parallelSpmvSmV1 $1 pgi
../runClusterTest.sh parallelSpmvSmV2 $1 pgi

cd ../

