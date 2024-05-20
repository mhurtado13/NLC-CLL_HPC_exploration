from SALib.analyze import sobol ##sobol is a type of SA analysis implemented in SALib
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Number of parameters
nparams = 2 

#Ranges for each parameter
#param_ranges = {
#    'cell_cell_repulsion_strength': [0, 75],
#    'cell_cell_adhesion_strength': [0, 2],
#    'relative_maximum_adhesion_distance': [0, 3.5],
#    'cell_BM_adhesion_strength': [0, 1],
#    'speed': [0, 1],
#    'migration_bias': [0, 1],
#    'secretion_rate': [0, 1],
#    'fluid_change_rate': [0, 100]
#}

param_ranges = {
    'cell_cell_repulsion_strength': [1, 10],
    'cell_cell_adhesion_strength': [0, 0.6]
}

#Define the problem for SALib
problem = {
    'num_vars': nparams,
    'names': list(param_ranges.keys()),
    'bounds': [[1, 10], [0, 0.6]] 
}

#Samples
samples_sobol = np.loadtxt('../data_output/sobol_samples.csv', delimiter=",", skiprows=1)

#Read output of simulation
output = np.loadtxt('../data_output/viability.csv', delimiter=",", skiprows=1)

#Take the median of output (as it's time series data)
output_median = np.median(output, axis=0)

# Sobol analysis
Si = sobol.analyze(problem, output_median, print_to_console=True, calc_second_order=True) 

print('Sobol Analysis Results:')
for i, param in enumerate(problem['names']):
   print(f"{param}: S1={Si['S1'][i]:.3f}, ST={Si['ST'][i]:.3f}")

# Plotting results
axes = Si.plot()
axes[0].set_yscale('log')
fig = plt.gcf()  # get current figure
fig.set_size_inches(10, 4)
plt.tight_layout()
