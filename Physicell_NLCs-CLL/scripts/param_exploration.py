import pandas as pd
import sys
import numpy as np

#Simulations = 20 variables x 5 values x 10 replicates = 1000 simulations

num_nodes = int(sys.argv[1])

input = {'uptake_rate_cancer': 1.0, 'secretion_rate_cancer':1.0, 'speed_cancer': 1.0, 'transformation_rate_cancer': 3e-5, 'relative_maximum_adhesion_distance_cancer': 0, 'cell_adhesion_affinity_cancer': 1,
         'speed_monocytes': 1.0, 'transformation_rate_monocytes': 95e-6, 'phagocytosis_rate_monocytes':25e-2, 'speed_macrophages': 1.0, 'phagocytosis_rate_macrophages':92e-2, 
         'attack_rate_macrophages': 5e-2, 'relative_maximum_adhesion_distance_macrophages': 0, 'cell_adhesion_affinity_macrophages': 1, 'speed_NLCs': 1.0, 'phagocytosis_rate_NLCs':4e-2, 
         'secretion_rate_apoptotic': 1.0, 'speed_apoptotic': 1.0, 'transformation_rate_apoptotic': 5e-05, 'secretion_rate_dead': 1.0}

default_values = list(input.values())

explore_values = list([1e-6, 1e-3, 1e-1, 2, 5])

def reset_values(data, values_def):        
    for i, key in enumerate(data.keys()):
        data[key] = values_def[i]

x = []

for parameter in input.keys():
    for i in explore_values:
        input[parameter] = i
        x.append(tuple(input.values()))
        reset_values(input, default_values)

#Subspaces for running in different nodes
rows = int(len(x)/num_nodes)

# Loop over the number of output files to generate
for i in range(num_nodes):
    # Calculate the start and end indices for the current output file
    start_idx = i * rows
    end_idx = (i + 1) * rows

    # If this is the last file, include any remaining rows
    if i == num_nodes - 1:
        end_idx = len(x)

    # Extract the rows for the current output file
    subset = x[start_idx:end_idx]

    thread_params = pd.DataFrame(subset)
    filename = f'data_output/Parameter_exploration/samples/Samples_{i}.csv'
    thread_params.to_csv(filename, index=False)





