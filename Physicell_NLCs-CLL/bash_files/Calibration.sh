#!/bin/bash
#SBATCH --job-name=NLC-CLL-calibration
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=120
#SBATCH --mem=128G
#SBATCH --partition=unlimitq
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=marcelo.hurtado@inserm.fr
#SBATCH -o /home/mhurtado/work/NLC-CLL_HPC_exploration/logs/Calibration/%x.o
#SBATCH -e /home/mhurtado/work/NLC-CLL_HPC_exploration/logs/Calibration/%x.e

# Define the number of tasks running in parallel
NUM_TASKS=2

# Define the number replicates for bootstrapping 
NUM_REPLICATES=2

# Define the population size for genetic algorithm (e.g. 500)
POP_SIZE=500

# Define the number of generations for genetic algorithm (e.g. 100)
NUM_GENERATION=100

#Specify which node to use
NODE=1

# Calibration of Physicell model using NSGA-II
python -u scripts/src/param_calibration.py $NUM_TASKS $NUM_REPLICATES $POP_SIZE $NUM_GENERATION $NODE
