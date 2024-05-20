import sys
import pandas as pd
from multiprocessing.pool import ThreadPool
from multiprocessing import Pool
from sampling import sobol_sampling
from model_simulation import run_model

nsamples = int(sys.argv[1])
num_tasks = int(sys.argv[2])
n_replicates = int(sys.argv[3]) 

pool = ThreadPool(num_tasks) 

#Load samples from Sobol
samples_sobol_all = sobol_sampling(nsamples)
samples_sobol_all = tuple(samples_sobol_all.itertuples(index=False, name=None))

thread_params = []

for i, param in enumerate(samples_sobol_all):
    thread_id = i % num_tasks + 1
    thread_params.append((thread_id,) + param)

params = [(("config/NLC_CLL.xml", n_replicates) + thread_params[contador]) for contador in range(len(thread_params))]
results = pool.starmap(run_model, params)

pool.close()
pool.join()

print("Pool closed")
print("Everything done! Results are saved in the ./data_output folder")

#Initialize viability and concentration vectors with first results
viability = results[0][0]
concentration = results[0][1]

for i in range(1, len(results)):
    via, conc = results[i]
    viability = pd.concat([viability, via], axis=1, ignore_index=True) #concatenating in the same order as values
    concentration = pd.concat([concentration, conc], axis=1, ignore_index=True) #concatenating in the same order as values

viability.to_csv('data_output/viability_SA.csv', index=False, header=True)
concentration.to_csv('data_output/concentration_SA.csv', index=False, header=True)