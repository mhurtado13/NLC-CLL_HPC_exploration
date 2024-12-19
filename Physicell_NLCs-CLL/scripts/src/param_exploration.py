import pandas as pd
import sys
import numpy as np

#Simulations = 33 variables x 8 values x 5 replicates = 1320 simulations

num_nodes = int(sys.argv[1])

input = {'cell_cell_repulsion_strength_cancer': 30, 'speed_cancer': 2,
         'secretion_rate_cytokines_cancer': 1e-1, 'uptake_rate_antiapoptotic_cancer': 1, 
         'transformation_rate_apoptotic_cancer': 8e-5, 'speed_monocytes': 1e-1,
         'uptake_rate_cytokines_monocytes': 1, 'uptake_rate_stress_monocytes': 1,
         'phagocytosis_rate_apoptotic_monocytes': 25e-4, 'phagocytosis_rate_dead_monocytes': 25e-4,
         'transformation_rate_macrophages_monocytes': 2e-12, 'transformation_rate_NLCs_monocytes': 3e-15,
         'speed_macrophages': 1e-1, 'uptake_rate_cytokines_macrophages': 1, 'uptake_rate_stress_macrophages': 1,
         'attack_rate_cancer_macrophages': 5e-2, 'damage_rate_macrophages': 10e-3, 'phagocytosis_rate_apoptotic_macrophages': 92e-4,
         'phagocytosis_rate_dead_macrophages': 92e-4, 'transformation_rate_NLCs_macrophages': 3e-16,
         'cell_cell_adhesion_strength_NLCs': 1, 'attachment_rate_NLCs': 0.001, 'detachment_rate_NLCs': 0.0001, "speed_NLCs": 1e-1,
         'secretion_rate_antiapoptotic_NLCs': 1, 'uptake_rate_stress_NLCs': 1, 'phagocytosis_rate_apoptotic_NLCs': 4e-4,
         'phagocytosis_rate_dead_NLCs': 4e-4, 'speed_apoptotic': 2, 'secretion_rate_cytokines_apoptotic': 1e-3,
         'secretion_rate_stress_apoptotic': 1e-1, 'transformation_rate_dead_apoptotic': 5e-5, 'secretion_rate_stress_dead': 1}

default_values = list(input.values())

explore_values = list([0, 1e-8, 1e-6, 1e-4, 1e-2, 2, 4, 6])

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





