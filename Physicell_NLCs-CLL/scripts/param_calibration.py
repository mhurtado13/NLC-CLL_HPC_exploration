import sys
import numpy as np
from multiprocessing import Pool
from pymoo.core.problem import Problem
from multiprocessing.pool import ThreadPool
from model_simulation import run_model

num_tasks = int(sys.argv[1]) 
n_replicates = int(sys.argv[2])
pop_size = int(sys.argv[3])
n_gen = int(sys.argv[4])

pool = ThreadPool(num_tasks) 

experimental = np.loadtxt('../Netlogo_NLCs-CLL/filtered_fused_9patients.csv', delimiter=",", skiprows=1)
viability_exp = experimental[:,1]
concentration_exp = experimental[:,2]

class calibrationProb(Problem):
    def __init__(self):
        super().__init__(n_var = 12, #12 parameters
                       n_obj = 2, #2 objective functions: viability and concentration
                       xl = np.array([0.9, 0.9, 4e-5, 0.9, 24e-2, 0.9, 91e-2, 0.9, 0.9, 3e-2, 1e-3, 0.9]), #lower bounds
                       xu = np.array([1.2, 1.2, 6e-5, 1.2, 26e-2, 1.2, 93e-2, 1.2, 1.2, 5e-2, 4e-3, 1.2])) #upper bounds
        
    def _evaluate(self, x, out):
        
        x = tuple(tuple(row) for row in x.tolist())
            
        thread_params = []

        for i, param in enumerate(x):
            thread_id = i % num_tasks + 1
            thread_params.append((thread_id,) + param)

        params = [(("config/NLC_CLL.xml", n_replicates) + thread_params[i]) for i in range(pop_size)]
        results = pool.starmap(run_model, params)

        #Objective functions
        obj1 = []
        obj2 = []
        for i in range(pop_size):
            viability, concentration = results[i]
            #RMSE of viability
            rmse_viability = np.sqrt(np.sum((viability - viability_exp)**2) / 10) #10 is the total of time points
            obj1.append(rmse_viability)
            #RMSE of concentration
            rmse_concentration = np.sqrt(np.sum((concentration - concentration_exp)**2) / 10) #10 is the total of time points
            obj2.append(rmse_concentration)

        #Stacking objectives to "F" 
        out["F"] = np.column_stack([obj1, obj2])


NLC_problem = calibrationProb()

from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.optimize import minimize
from pymoo.termination import get_termination

algorithm_nsga = NSGA2(pop_size=pop_size)

termination = get_termination("n_gen", n_gen)

res = minimize(NLC_problem,
               algorithm_nsga,
               termination,
               seed=1,
               verbose=True)

pool.close()
pool.join()

print(res.X)
print(res.F)

np.savetxt('data_output/Space_values.csv', res.X, delimiter=",")
np.savetxt('data_output/Objective_values.csv', res.F, delimiter=",")
