#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 01:30:00
#SBATCH --constraint=avx
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=e.maassen@uvt.nl

module load stopos
module load sara-batch-resources

# put the name of the pool in the environment variable
# STOPOS_POOL

export STOPOS_POOL= pool


# determine the number of cores available on the node the job is running on:
ncores=`sara-get-num-cores`

# -----------------------------------
# Setting up and running the simulation


for ((i=1; i<=ncores; i++)) ; do
(
	# Getting the next parameters from the pool
	stopos next

	# Checking if the parameters pool is empty
	if [ "$STOPOS_RC" != "OK" ]; then
		break
	fi

	pos="$STOPOS_VALUE"
    
	# Removing the used parameter from the pool
	stopos remove

	echo
	echo "Running the simulation for ${pos} position..."
	Rscript idealworld-lisa.R "${pos}"

) &
done
wait