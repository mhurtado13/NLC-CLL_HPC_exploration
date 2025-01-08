from SALib.sample import sobol
import pandas as pd
import numpy as np

def sobol_sampling(nsamples):
    
    
    ##number of inputs = 14

    default_values = {'cell_cell_repulsion_strength_cancer': 30, 'speed_cancer': 2, 'phagocytosis_rate_dead_monocytes': 25e-4,
                      'transformation_rate_NLCs_monocytes': 3e-15, 'uptake_rate_cytokines_macrophages': 1,
                      'attack_rate_cancer_macrophages': 5e-2, 'phagocytosis_rate_apoptotic_macrophages': 92e-4,
                      'phagocytosis_rate_dead_macrophages': 92e-4, 'transformation_rate_NLCs_macrophages': 3e-16,
                      'cell_cell_adhesion_strength_NLCs': 1, 'attachment_rate_NLCs': 0.001, 
                      'phagocytosis_rate_dead_NLCs': 4e-4, 'speed_apoptotic': 2, 'transformation_rate_dead_apoptotic': 5e-5}
    
    
    #Define the problem for SALib
    problem = {
        'num_vars': len(default_values),
        'names': default_values.keys(),
        'bounds': np.array([[10, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                            [60, 4, 25e-2, 6e-15, 2, 10e-2, 92e-2, 92e-2, 6e-16, 2, 0.002, 4e-2, 3, 10e-5]]).T 
    }

    ## number of inputs = 128

    #Generate sobol samples: N*(2D + 2) samples where N: nsamples and D: number of inputs. For sobol analysis N needs to be a power of 2, if not it will give error 
    sobol_samples = sobol.sample(problem, nsamples) 

    param_names = list(default_values.keys())
    sobol_samples = pd.DataFrame(sobol_samples, columns=param_names)

    return sobol_samples


