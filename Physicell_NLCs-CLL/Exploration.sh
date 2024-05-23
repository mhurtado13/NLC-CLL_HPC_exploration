#!/bin/bash
#SBATCH --job-name=NLC-CLL-exploration
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=100
#SBATCH --mem=256G
#SBATCH --partition=workq
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=marcelo.hurtado@inserm.fr
#SBATCH -o /home/mhurtado/work/NLC-CLL_HPC_exploration/logs/Exploration/%x.o
#SBATCH -e /home/mhurtado/work/NLC-CLL_HPC_exploration/logs/Exploration/%x.e

# Parameter exploration
#python scripts/param_exploration.py $NUM_NODES $NUM_TASKS 

#Define file of samples to run exploration
FILE=data_output/Parameter_exploration/Samples_1.csv
#Define the number of tasks running in parallel
NUM_TASKS=8
#Define the number replicates for bootstrapping 
NUM_REPLICATES=10
#Specify which node to use
NODE=1

python scripts/run_parallel.py $FILE $NUM_TASKS $NUM_REPLICATES $NODE

