#!/bin/bash
#SBATCH --job-name=NLC-CLL-calibration
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=120
#SBATCH --mem=128G
#SBATCH --partition=unlimitq
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=marcelo.hurtado@inserm.fr
#SBATCH -o /home/mhurtado/work/ABM_NLC-CLL_cluster/logs/Calibration/%x.o
#SBATCH -e /home/mhurtado/work/ABM_NLC-CLL_cluster/logs/Calibration/%x.e

# Define the number of tasks running in parallel
NUM_TASKS=120

# Define the number replicates for bootstrapping 
NUM_REPLICATES=10

# Define the population size for genetic algorithm
POP_SIZE=500

# Define the number of generations for genetic algorithm
NUM_GENERATION=100

# Calibration of Physicell model using NSGA-II
python scripts/param_calibration.py $NUM_TASKS $NUM_REPLICATES $POP_SIZE $NUM_GENERATION
