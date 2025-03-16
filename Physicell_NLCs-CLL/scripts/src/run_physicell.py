from model_simulation import simulate_model
import numpy as np
from pymoo.util.nds.non_dominated_sorting import NonDominatedSorting
import matplotlib.pyplot as plt

# Load the objective values
objective_values = np.loadtxt('data_output/Calibration/Objective_values.csv', delimiter=",")

# Extract objectives
obj1 = objective_values[:, 0]  # RMSE Viability
obj2 = objective_values[:, 1]  # RMSE Concentration

# Identify Pareto front
nds = NonDominatedSorting().do(objective_values, only_non_dominated_front=True)

pareto_obj1 = obj1[nds]
pareto_obj2 = obj2[nds]

# Find best viability and concentration points
best_viability_idx = np.argmin(pareto_obj1)
best_concentration_idx = np.argmin(pareto_obj2)

# Load results
parameter_values = np.loadtxt('data_output/Calibration/Space_values.csv', delimiter=",")

# Get the corresponding parameter sets
pareto_params_viability = parameter_values[best_viability_idx]
pareto_params_concentration = parameter_values[best_concentration_idx]

# Print and save results
print("Non-Dominated Parameter Sets:")
print("Viability:")
print(pareto_params_viability)
print("Concentration:")
print(pareto_params_concentration)

params = ("config/NLC_CLL.xml", 5, 1, 1) + tuple(pareto_params_viability)
results = simulate_model(*params)
viability, concentration,_ = results

viability_name = f'data_output/Calibration/viability_calibrated.csv'
concentration_name = f'data_output/Calibration/concentration_calibrated.csv'

viability.to_csv(viability_name, index=False, header=True)
concentration.to_csv(concentration_name, index=False, header=True)

print("Viability:")
print(viability)
print("Concentration:")
print(concentration)