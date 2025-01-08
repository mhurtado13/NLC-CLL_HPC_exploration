#!/bin/bash
#SBATCH --job-name=NLC-CLL-exploration-node2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --partition=unlimitq
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=marcelo.hurtado@inserm.fr
#SBATCH -o /home/mhurtado/work/NLC-CLL_HPC_exploration/Physicell_NLCs-CLL/logs/Exploration/%x.o
#SBATCH -e /home/mhurtado/work/NLC-CLL_HPC_exploration/Physicell_NLCs-CLL/logs/Exploration/%x.e

#Define file of samples to run exploration
FILE="data_output/Parameter_exploration/samples/Samples_1.csv"
#Define the number of tasks running in parallel
NUM_TASKS=32
#Define the number replicates for bootstrapping 
NUM_REPLICATES=5
#Specify which node to use
NODE=2

python -u scripts/src/run_model.py $FILE $NUM_TASKS $NUM_REPLICATES $NODE #flag -u ensures printing log messages for real-time debugging

#Before running this, run in the console: python scripts/param_exploration.py $NUM_NODES
