#!/bin/bash

#======================================================
#
# Job script for running ROMS on multiple nodes
#
#======================================================

#======================================================
# Propogate environment variables to the compute node
#SBATCH --export=ALL
#
# Run in the standard partition (queue)
#SBATCH --partition=standard
#
# Specify project account
#SBATCH --account=chen-dse
#
# No. of tasks required
#SBATCH --ntasks=5 --nodes=1
#
# Distribute processes in round-robin fashion for load balancing
#SBATCH --distribution=cyclic
#
# Specify (hard) runtime (HH:MM:SS)
#SBATCH --time=120:00:00
#
# Job name
#SBATCH --job-name=IBM
#
# Output file
#SBATCH --output=IBM-%j.out
#======================================================

module purge
module load netcdf/intel-2020.4/4.9.0
module load netcdf-fortran/intel-2020.4/4.5.4

#export LD_LIBRARY_PATH=/opt/software/netcdf/intel-2018.2/4.6.1/lib:$LD_LIBRARY_PATH
#======================================================
# Prologue script to record job details
# Do not change the line below
#======================================================
/opt/software/scripts/job_prologue.sh  
#------------------------------------------------------

mpirun -np $SLURM_NPROCS ./IBM > out

#======================================================
# Epilogue script to record job endtime and runtime
# Do not change the line below
#======================================================
/opt/software/scripts/job_epilogue.sh 
#------------------------------------------------------
