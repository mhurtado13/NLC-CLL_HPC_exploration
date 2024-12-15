import sys
import pandas as pd
from sampling import sobol_sampling

nsamples = int(sys.argv[1])
num_nodes = int(sys.argv[2])

#Load samples from Sobol
samples_sobol_all = sobol_sampling(nsamples)

#N*(2D + 2) N = nsamples D = inputs 

#Subspaces for running in different nodes
rows = int(len(samples_sobol_all)/num_nodes)

# Loop over the number of output files to generate
for i in range(num_nodes):
    # Calculate the start and end indices for the current output file
    start_idx = i * rows
    end_idx = (i + 1) * rows

    # If this is the last file, include any remaining rows
    if i == num_nodes - 1:
        end_idx = len(samples_sobol_all)

    # Extract the rows for the current output file
    subset = samples_sobol_all[start_idx:end_idx]

    thread_params = pd.DataFrame(subset)
    filename = f'data_output/Sensitivity_analysis/samples/Samples_{i}.csv'
    thread_params.to_csv(filename, index=False)