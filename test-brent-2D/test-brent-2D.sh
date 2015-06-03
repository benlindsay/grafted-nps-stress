#!/bin/sh
#PBS -l nodes=1:ppn=3
#PBS -l walltime=148:00:00,mem=2gb

cd $PBS_O_WORKDIR

while read machine
do
  echo $machine
done < $PBS_NODEFILE > nodes

NO_OF_NODES=`cat $PBS_NODEFILE | egrep -v '^#'\|'^$' | wc -l | awk '{print $1}'`
NODE_LIST=`cat $PBS_NODEFILE`

module load intel-11.1
module load openmpi-1.4.3-intel

time mpirun -np $NO_OF_NODES -machinefile nodes ../a.out > LOG
