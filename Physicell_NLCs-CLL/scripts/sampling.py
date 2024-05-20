from SALib.sample import sobol
import pandas as pd
import numpy as np

def sobol_sampling(nsamples):

    default_values = {'uptake_rate_cancer': 1.0, 'speed_cancer': 1.0, 'transformation_rate_cancer': 5e-5,
                  'speed_monocytes':1.0, 'dead_phagocytosis_rate_monocytes':25e-2, 'speed_macrophages':1.0,
                  'dead_phagocytosis_rate_macrophages':92e-2, 'secretion_rate_NLCs':1.0, 'speed_NLCs':1.0,
                  'dead_phagocytosis_rate_NLCs':4e-2, 'death_rate_apoptotic':3e-3, 'secretion_rate_apoptotic': 1.0}
    ## number of inputs = 12

    #Define the problem for SALib
    problem = {
        'num_vars': len(default_values),
        'names': default_values.keys(),
        'bounds': np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [2.0, 2.0, 10e-5, 2.0, 50e-2, 2.0, 184e-2, 2.0, 2.0, 8e-2, 6e-3, 2]]).T # 100% percent of variation from default values 
    }

    #Generate sobol samples: N*(2D + 2) samples where N: nsamples and D: number of inputs 
    sobol_samples = sobol.sample(problem, nsamples) 

    param_names = list(default_values.keys())
    sobol_samples = pd.DataFrame(sobol_samples, columns=param_names)

    return sobol_samples

    #Save output
    #sobol_samples.to_csv('data_output/sobol_samples.csv', index=False)

