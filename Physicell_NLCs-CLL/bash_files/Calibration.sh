#!/bin/bash
#SBATCH --job-name=NLC-CLL-calibration
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --partition=long
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=marcelo.hurtado@inserm.fr
#SBATCH -o ~/NLC-CLL_HPC_exploration/Physicell_NLCs-CLL/logs/Calibration/%x.o
#SBATCH -e ~/NLC-CLL_HPC_exploration/Physicell_NLCs-CLL/logs/Calibration/%x.e

# Define the number of tasks running in parallel
NUM_TASKS=32

# Define the number replicates for bootstrapping 
NUM_REPLICATES=1

# Define the population size for genetic algorithm (e.g. 500)
POP_SIZE=200

# Define the number of generations for genetic algorithm (e.g. 100)
NUM_GENERATION=100

#Specify which node to use
NODE=1

# Calibration of Physicell model using NSGA-II
python -u scripts/src/param_calibration.py $NUM_TASKS $NUM_REPLICATES $POP_SIZE $NUM_GENERATION $NODE
