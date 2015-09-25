#!/bin/bash

#SBATCH -J nobrent4.601-bcheck     # job name
#SBATCH -o nobrent4.601-bcheck.o%j # output and error file name (%j expands to jobID)
#SBATCH -n 8               # total number of mpi tasks requested
#SBATCH -p normal          # queue (partition) -- normal, development, etc.
#SBATCH -t 48:00:00        # run time (hh:mm:ss)
#SBATCH --mail-user=benjlindsay@gmail.com
#SBATCH --mail-type=begin  # email me when the job starts
#SBATCH --mail-type=end    # email me when the job finishes

time ibrun ../debugbrent-1.1-2.exe > LOG
