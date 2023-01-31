#!/bin/sh

#SBATCH --job-name=solution_05
#SBATCH --partition=nu-wide
#SBATCH --ntasks=4
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=64
#SBATCH --time=03:00:00

#SBATCH --output=output_%j.out
#SBATCH --error=output_%j.err

PROGRAM=solution_05


mpirun -n 4 ./${PROGRAM}
