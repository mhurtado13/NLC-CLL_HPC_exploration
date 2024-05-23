#!/bin/bash
#SBATCH --job-name=NLC-CLL-exploration
#SBATCH --nodes=2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=100
#SBATCH --mem=256G
#SBATCH --partition=workq
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=marcelo.hurtado@inserm.fr
#SBATCH -o /home/mhurtado/work/NLC-CLL_HPC_exploration/logs/Exploration/%x.o
#SBATCH -e /home/mhurtado/work/NLC-CLL_HPC_exploration/logs/Exploration/%x.e

# Define the number of tasks running in parallel
NUM_TASKS=100

# Define the number replicates for bootstrapping 
NUM_REPLICATES=10

# Parameter exploration
python scripts/param_exploration.py $NUM_TASKS $NUM_REPLICATES 

#Simulations = 12 variables x 20 values x 10 replicates = 2400 simulations
