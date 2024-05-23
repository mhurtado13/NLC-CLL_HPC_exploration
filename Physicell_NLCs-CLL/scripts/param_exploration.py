import pandas as pd
import sys
import numpy as np

#Simulations = 12 variables x 20 values x 10 replicates = 2400 simulations

num_nodes = int(sys.argv[1])
num_tasks = int(sys.argv[2])

input = {'uptake_rate_cancer': 1.0, 'speed_cancer': 1.0, 'transformation_rate_cancer': 5e-5,
                  'speed_monocytes':1.0, 'dead_phagocytosis_rate_monocytes':25e-2, 'speed_macrophages':1.0,
                  'dead_phagocytosis_rate_macrophages':92e-2, 'secretion_rate_NLCs':1.0, 'speed_NLCs':1.0,
                  'dead_phagocytosis_rate_NLCs':4e-2, 'death_rate_apoptotic':3e-3, 'secretion_rate_apoptotic': 1.0}

#death rate apoptotic cells
#number of initial apoptotic cells 
#number of initial CLL cells

default_values = list(input.values())

explore_values = list(np.round(np.linspace(0, 1, 20), 2))
#explore_values = [0,1]

def reset_values(data, values_def):        
    for i, key in enumerate(data.keys()):
        data[key] = values_def[i]

x = []

for parameter in input.keys():
    for i in explore_values:
        input[parameter] = i
        x.append(tuple(input.values()))
        reset_values(input, default_values)

# 12 * 20 = 240 parameters

#Subspaces for running in different nodes
rows = int(len(x)/num_nodes)

# Loop over the number of output files to generate
for i in range(num_nodes):
    thread_params = []
    # Calculate the start and end indices for the current output file
    start_idx = i * rows
    end_idx = (i + 1) * rows

    # If this is the last file, include any remaining rows
    if i == num_nodes - 1:
        end_idx = len(x)

    # Extract the rows for the current output file
    subset = x[start_idx:end_idx]

    for k, param in enumerate(subset):
        thread_id = k % num_tasks + 1
        thread_params.append((thread_id,) + param)

    thread_params = pd.DataFrame(thread_params)
    filename = f'data_output/Parameter_exploration/Samples_{i}.csv'
    thread_params.to_csv(filename, index=False)





