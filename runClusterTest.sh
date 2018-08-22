#!/bin/bash
if [ "$#" -ne 3 ]
then
  echo "Usage: $0 program matrixSize compiler"
  exit 1
fi

ppn="$(($2 / 4))"


tempFilename='anyTempFileNameWillWork.txt'
outputFilename=$1.txt

nloops=5

# Determining MPI implementation and binding options #
MPI=`mpiexec --version | head -1 | awk '{print $1}' `

if [ "$MPI" == "HYDRA" ]; then
    echo "MPICH"
    bindings="--bind-to socket"
    export HYDRA_TOPO_DEBUG=1
    export MPIR_CVAR_CH3_PORT_RANGE=10000:10100
elif [ "$MPI" == "Intel(R)" ]; then
    echo "Intel MPI"
    bindings="-genv I_MPI_PIN_DOMAIN=core -genv I_MPI_PIN_ORDER=spread -genv I_MPI_DEBUG=4 -genv I_MPI_FABRICS=shm:tcp"
elif [ "$MPI" == "mpiexec" ]; then
    echo "open-mpi"
    bindings="--bind-to core --report-bindings --mca btl_tcp_if_exclude docker0,127.0.0.1/8"
fi
# end of Determining MPI implementation and binding options #

npt=`grep -c ^processor /proc/cpuinfo`
numaNodes=`lscpu | grep "NUMA node(s):" | awk '{}{print $3}{}'`
tpc=`lscpu | grep "Thread(s) per core:" | awk '{}{print $4}{}'`
np="$(($npt / $tpc))"
npps="$(($np / $numaNodes))"
npm1="$(($np - 1))"


if [ -n "$PGI" ]; then
    echo "Pgi Compiler"
elif [ -n "$INTEL_LICENSE_FILE" ]; then
    echo "Intel Compiler"
    #np=15
    #npps="$(($np / $numaNodes))"
    #npm1="$(($np - 1))"
else
    echo "Gnu Compiler"
fi

rm -f $tempFilename

for j in  `seq 1 $nloops`; do
    echo run number: $j
    if [ "$MPI" == "Intel(R)" ]; then
        mpiexec $bindings -hosts stout,koelsch,dunkel,porter  -n $2  -ppn $ppn $1 | grep taken >>  $tempFilename
    elif [ "$MPI" == "mpiexec" ]; then
        mpiexec $bindings -host stout:$ppn,koelsch:$ppn,dunkel:$ppn,porter:$ppn -n $2 $1 | grep taken >>  $tempFilename
    fi
done

mkdir -p ../plots/cluster/$3/$2

cat $tempFilename | awk 'BEGIN{} {print $5} END{}' | sort -n  | head -1 > ../plots/cluster/$3/$2/$outputFilename

rm $tempFilename

