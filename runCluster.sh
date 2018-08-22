#!/bin/bash
if [ "$#" -ne 1 ]
then
  echo "Usage: $0 matrixSize"
  exit 1
fi

cd buildGnu

../runClusterTest.sh test_exxon_readerNormal $1 gnu
../runClusterTest.sh test_exxon_readerSm     $1 gnu
#../runClusterTest.sh test_exxon_readerBSR    $1 gnu

cd ../

source setIcc intel64
source setImpi

cd buildIntel

../runClusterTest.sh test_exxon_readerNormal $1 intel
../runClusterTest.sh test_exxon_readerSm     $1 intel
#../runClusterTest.sh test_exxon_readerBSR    $1 intel


cd ../

source setPgi 18.x
source setPgiMpi 18.x

cd buildPgi

../runClusterTest.sh test_exxon_readerNormal $1 pgi
../runClusterTest.sh test_exxon_readerSm     $1 pgi
#../runClusterTest.sh test_exxon_readerBSR    $1 pgi

cd ../

