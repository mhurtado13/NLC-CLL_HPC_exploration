from multiprocessing import Pool
from model_simulation import simulate_model
from multiprocessing.pool import ThreadPool
import pandas as pd
import sys
import shutil
import os

file_params = str(sys.argv[1])
num_tasks = int(sys.argv[2])
n_replicates = int(sys.argv[3])
n_node = int(sys.argv[4])

pool = ThreadPool(num_tasks) 
results = []

thread_params = pd.read_csv(file_params)
thread_params = [tuple(row) for row in thread_params.itertuples(index=False, name=None)]
params = [(("config/NLC_CLL.xml", n_replicates, n_node, contador+1) + thread_params[contador]) for contador in range(len(thread_params))]
res = pool.starmap(simulate_model, params)
results.extend(res)

pool.close()
pool.join()

print("Pool closed")
print("Everything done in node " + str(n_node) + "! Results are saved in the ./data_output folder")

#Clean memory
shutil.rmtree("output/output_node_" + str(n_node)) #Remove output folder for each node when it is done 
shutil.rmtree("config/config_node_" + str(n_node)) #Remove config folder for each node when it is done

#Initialize viability and concentration vectors with first results
viability = results[0][0]
concentration = results[0][1]
param_error = []
for i in range(1, len(results)):
    via, conc, err = results[i]
    viability = pd.concat([viability, via], axis=1, ignore_index=True) #concatenating in the same order as explore_values
    concentration = pd.concat([concentration, conc], axis=1, ignore_index=True) #concatenating in the same order as explore_values
    param_error.extend(err)

viability_name = f'data_output/Sensitivity_analysis/results/viability_{n_node}.csv'
concentration_name = f'data_output/Sensitivity_analysis/results/concentration_{n_node}.csv'
errors_name = f'data_output/Sensitivity_analysis/results/Physicell_errors_{n_node}.csv'

viability.to_csv(viability_name, index=False, header=True)
concentration.to_csv(concentration_name, index=False, header=True)

if param_error:
    param_error_df = pd.DataFrame(param_error)
    param_error_df.to_csv(errors_name, index=False, header=True)
    print("Physicell errors in node " + str(n_node) + " for parameters:\n" + str(param_error))
else:
    print("All parameters were evaluated succesfully in node " + str(n_node))