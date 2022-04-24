#!/bin/bash

#SBATCH --job-name=HybridModel
#SBATCH -p general
#SBATCH -o Hybrid_%j.txt
#SBATCH -e Hybrid_%j.err
#SBATCH --ntasks=47
#SBATCH --time=0-1:00:00
#SBATCH --mem=16G

module load python/3.6.11
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

srun --cpu-bind=sockets python Single_Run.py