from SALib.sample import sobol
import pandas as pd
import numpy as np

def sobol_sampling(nsamples):
    
    default_values = {'uptake_rate_cancer': 1.0, 'secretion_rate_cancer':1.0, 'speed_cancer': 1.0, 'transformation_rate_cancer': 3e-5, 'relative_maximum_adhesion_distance_cancer': 0, 'cell_adhesion_affinity_cancer': 1,
                      'speed_monocytes': 1.0, 'transformation_rate_monocytes': 95e-6, 'phagocytosis_rate_monocytes':25e-2, 'speed_macrophages': 1.0, 'phagocytosis_rate_macrophages':92e-2, 
                      'attack_rate_macrophages': 5e-2, 'relative_maximum_adhesion_distance_macrophages': 0, 'cell_adhesion_affinity_macrophages': 1, 'speed_NLCs': 1.0, 'phagocytosis_rate_NLCs':4e-2, 
                      'secretion_rate_apoptotic': 1.0, 'speed_apoptotic': 1.0, 'transformation_rate_apoptotic': 5e-05, 'secretion_rate_dead': 1.0}

    ## number of inputs = 20

    #Define the problem for SALib
    problem = {
        'num_vars': len(default_values),
        'names': default_values.keys(),
        'bounds': np.array([[0, 0, 0, 3e-7, 0, 0, 0, 95e-7, 25e-4, 0, 92e-4, 5e-7, 0, 0, 0, 4e-4, 0, 0, 5e-4, 0],
                            [3, 3, 3, 3e-3, 3, 3, 3, 95e-3, 25e-1, 3, 92e-1, 5e-3, 3, 3, 3, 4e-1, 3, 3, 5e-1, 3]]).T 
    }

    #Generate sobol samples: N*(2D + 2) samples where N: nsamples and D: number of inputs 
    sobol_samples = sobol.sample(problem, nsamples) 

    param_names = list(default_values.keys())
    sobol_samples = pd.DataFrame(sobol_samples, columns=param_names)

    return sobol_samples

    #Save output
    #sobol_samples.to_csv('data_output/sobol_samples.csv', index=False)

