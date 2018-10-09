#!/bin/bash
#   Crearte soft links with the programs into the same directory 
#   Do this for each compiler
#
if [ "$#" -ne 1 ]
then
  echo "Usage: $0 matrixName"
  exit 1
fi
matrix=$1

cd buildGnu

../runClusterTest.sh parallelSpmv   $matrix gnu
../runClusterTest.sh parallelSpmvSm $matrix gnu

cd ../

source setIcc intel64
source setImpi

cd buildIntel

../runClusterTest.sh parallelSpmv   $matrix intel
../runClusterTest.sh parallelSpmvSm $matrix intel


cd ../

source setPgi 18.x
source setPgiMpi 18.x

cd buildPgi

../runClusterTest.sh parallelSpmv   $matrix pgi
../runClusterTest.sh parallelSpmvSm $matrix pgi

cd ../

