#!/bin/bash
#SBATCH --job-name=NLC-CLL-exploration
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=120
#SBATCH --mem=128G
#SBATCH --partition=workq
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=marcelo.hurtado@inserm.fr
#SBATCH -o /home/mhurtado/work/ABM_NLC-CLL_cluster/logs/Exploration/%x.o
#SBATCH -e /home/mhurtado/work/ABM_NLC-CLL_cluster/logs/Exploration/%x.e

# Define the number of tasks running in parallel
NUM_TASKS=120

# Define the number replicates for bootstrapping 
NUM_REPLICATES=10

# Parameter exploration
python scripts/param_exploration.py $NUM_TASKS $NUM_REPLICATES 
