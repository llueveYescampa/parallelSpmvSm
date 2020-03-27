#!/bin/bash
if [ "$#" -lt 1 ]; then
  echo "Usage: $0 matrixName [numOfProcesses]"
  exit 1
fi
n=4
if [ "$#" -eq 2 ]; then
    n=$2
fi
echo "number of threads or processes: $n"
matrix=$1

nloops=5

tempFilename=$(hostname)'_anyTempFileNameWillWork.txt'
mkdir -p ../plots/$(hostname)'_'$n/Gnu/$matrix/
mkdir -p ../plots/$(hostname)'_'$n/Intel/$matrix/
mkdir -p ../plots/$(hostname)'_'$n/Pgi/$matrix/
#outputFilename=$1.txt

cd buildGnu

rm -f $tempFilename
for j in  `seq 1 $nloops`; do
    run raptor_Normal $matrix $n | grep by >>  $tempFilename
done
cat $tempFilename  | grep process | awk 'BEGIN{}   { printf("%f  %f\n", $7,$10)}  END{}' |  sort  -k1,1n -k2,2n |   head -1  > ../../plots/$(hostname)'_'$n/Gnu/$matrix/raptor_Normal.txt


rm -f $tempFilename
for j in  `seq 1 $nloops`; do
    run raptor_TAP $matrix $n | grep by >>  $tempFilename
done
cat $tempFilename  | grep process | awk 'BEGIN{}   { printf("%f  %f\n", $7,$10)}  END{}' |  sort  -k1,1n -k2,2n |   head -1  > ../../plots/$(hostname)'_'$n/Gnu/$matrix/raptor_TAP.txt


rm -f $tempFilename
for j in  `seq 1 $nloops`; do
    run parallelSpmv $matrix $n | grep by >>  $tempFilename
done
cat $tempFilename  | grep process | awk 'BEGIN{}   { printf("%f  %f\n", $7,$10)}  END{}' |  sort  -k1,1n -k2,2n |   head -1  > ../../plots/$(hostname)'_'$n/Gnu/$matrix/parallelSpmv.txt

rm -f $tempFilename
for j in  `seq 1 $nloops`; do
    run parallelSpmvSm $matrix $n | grep by >>  $tempFilename
done
cat $tempFilename  | grep process | awk 'BEGIN{}   { printf("%f  %f\n", $7,$10)}  END{}' |  sort  -k1,1n -k2,2n |   head -1  > ../../plots/$(hostname)'_'$n/Gnu/$matrix/parallelSpmvSm.txt

rm $tempFilename

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

cd ../buildIntel

rm -f $tempFilename
for j in  `seq 1 $nloops`; do
    run raptor_Normal $matrix $n | grep by >>  $tempFilename
done
cat $tempFilename  | grep process | awk 'BEGIN{}   { printf("%f  %f\n", $7,$10)}  END{}' |  sort  -k1,1n -k2,2n |   head -1  > ../../plots/$(hostname)'_'$n/Intel/$matrix/raptor_Normal.txt


rm -f $tempFilename
for j in  `seq 1 $nloops`; do
    run raptor_TAP $matrix $n | grep by >>  $tempFilename
done
cat $tempFilename  | grep process | awk 'BEGIN{}   { printf("%f  %f\n", $7,$10)}  END{}' |  sort  -k1,1n -k2,2n |   head -1  > ../../plots/$(hostname)'_'$n/Intel/$matrix/raptor_TAP.txt


rm -f $tempFilename
for j in  `seq 1 $nloops`; do
    run parallelSpmv $matrix $n | grep by >>  $tempFilename
done
cat $tempFilename  | grep process | awk 'BEGIN{}   { printf("%f  %f\n", $7,$10)}  END{}' |  sort  -k1,1n -k2,2n |   head -1  > ../../plots/$(hostname)'_'$n/Intel/$matrix/parallelSpmv.txt

rm -f $tempFilename
for j in  `seq 1 $nloops`; do
    run parallelSpmvSm $matrix $n | grep by >>  $tempFilename
done
cat $tempFilename  | grep process | awk 'BEGIN{}   { printf("%f  %f\n", $7,$10)}  END{}' |  sort  -k1,1n -k2,2n |   head -1  > ../../plots/$(hostname)'_'$n/Intel/$matrix/parallelSpmvSm.txt


rm -f $tempFilename
for j in  `seq 1 $nloops`; do
    run mklSpmv $matrix $n | grep by >>  $tempFilename
done
cat $tempFilename  | grep threads | awk 'BEGIN{}   { printf("%f  %f\n", $7,$10)}  END{}' |  sort  -k1,1n -k2,2n |   head -1  > ../../plots/$(hostname)'_'$n/Intel/$matrix/mklSpmv.txt

rm -f $tempFilename
for j in  `seq 1 $nloops`; do
    run parallelSpmvSmAlphaBeta $matrix $n | grep by >>  $tempFilename
done
cat $tempFilename  | grep process | awk 'BEGIN{}   { printf("%f  %f\n", $7,$10)}  END{}' |  sort  -k1,1n -k2,2n |   head -1  > ../../plots/$(hostname)'_'$n/Intel/$matrix/parallelSpmvSmAlphaBeta.txt



rm $tempFilename

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
cd ../buildPgi

rm -f $tempFilename
for j in  `seq 1 $nloops`; do
    run raptor_Normal $matrix $n | grep by >>  $tempFilename
done
cat $tempFilename  | grep process | awk 'BEGIN{}   { printf("%f  %f\n", $7,$10)}  END{}' |  sort  -k1,1n -k2,2n |   head -1  > ../../plots/$(hostname)'_'$n/Pgi/$matrix/raptor_Normal.txt


rm -f $tempFilename
for j in  `seq 1 $nloops`; do
    run raptor_TAP $matrix $n | grep by >>  $tempFilename
done
cat $tempFilename  | grep process | awk 'BEGIN{}   { printf("%f  %f\n", $7,$10)}  END{}' |  sort  -k1,1n -k2,2n |   head -1  > ../../plots/$(hostname)'_'$n/Pgi/$matrix/raptor_TAP.txt


rm -f $tempFilename
for j in  `seq 1 $nloops`; do
    run parallelSpmv $matrix $n | grep by >>  $tempFilename
done
cat $tempFilename  | grep process | awk 'BEGIN{}   { printf("%f  %f\n", $7,$10)}  END{}' |  sort  -k1,1n -k2,2n |   head -1  > ../../plots/$(hostname)'_'$n/Pgi/$matrix/parallelSpmv.txt

rm -f $tempFilename
for j in  `seq 1 $nloops`; do
    run parallelSpmvSm $matrix $n | grep by >>  $tempFilename
done
cat $tempFilename  | grep process | awk 'BEGIN{}   { printf("%f  %f\n", $7,$10)}  END{}' |  sort  -k1,1n -k2,2n |   head -1  > ../../plots/$(hostname)'_'$n/Pgi/$matrix/parallelSpmvSm.txt

rm $tempFilename
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



cd
echo "bye ..."

