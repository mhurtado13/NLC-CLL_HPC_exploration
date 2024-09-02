import numpy as np
import pandas as pd
import pcdl
import shutil
import os

#config_file = 'config/NLC_CLL.xml'
def collect(dir_output, config_file, node):
    #Create error folder
    error_folder = "error/error_node_" + str(node)
    os.makedirs(error_folder, exist_ok=True)

    mcdsts = pcdl.TimeSeries(dir_output, settingxml=config_file, verbose = False) 
    timesteps = mcdsts.get_mcds_list()

    #Extract positions corresponding to days 1-13
    positions = []
    for days in range(0,14):
        hours = 24*days
        positions.append(hours)

    #Initial CLL cells
    initial = timesteps[0].get_cell_df(states=1)
    CLL_initial = len(initial[(initial['cell_type']=="cancer")])
    apoptotic_initial = len(initial[(initial['cell_type']=="apoptotic")])
    dead_initial = len(initial[(initial['cell_type']=="dead")])

    alive = [CLL_initial]
    dead = [dead_initial]
    apoptotic = [apoptotic_initial]
    for i in range(1, len(positions)):
        step = timesteps[positions[i]].get_cell_df(states=1)
        number_alive = len(step[(step['cell_type']=='cancer')&(step['dead']==False)]) #step['dead'] is only a formality cause all cells are considered 'alive', 'dead' is another celltype for this model
        number_apoptotic = len(step[(step['cell_type']=='apoptotic')&(step['dead']==False)])
        number_dead = len(step[(step['cell_type']=='dead')&(step['dead']==False)])
        alive.append(number_alive)
        dead.append(number_dead)
        apoptotic.append(number_apoptotic)

    CLL_alive = pd.Series(alive, name="Cells_alive")
    CLL_apoptotic = pd.Series(apoptotic, name = "Cells_apoptotic")
    CLL_dead = pd.Series(dead, name = "Cells_dead")

    #viability at time t =  CLL alive at time t / (CLL alive + CLL apoptotic + CLL dead) at time t
    viability = []
    for i in range(len(CLL_alive)):
        try:
            # Calculate the number and append to viability list
            number = (CLL_alive[i] / (CLL_alive[i] + CLL_apoptotic[i] + CLL_dead[i])) * 100
            viability.append(number)
        except Exception as e:
            destination_file = os.path.join(error_folder, os.path.basename(config_file))
            shutil.copy(config_file, destination_file)
            number = 0  # Continue
            viability.append(number)

    ####Remove day 4, 5, 11, 12 because of experimental
    viability = np.delete(viability, [4,5,11,12], axis=0)

    viability = pd.Series(viability, name = "CLL viability")

    #concentration at time t =  CLL alive at time t / (CLL initial)
    concentration = []
    for i in range(len(CLL_alive)):
        try:
            number = (CLL_alive[i]/CLL_initial)*100
            concentration.append(number)
        except Exception as e:
            destination_file = os.path.join(error_folder, os.path.basename(config_file))
            shutil.copy(config_file, destination_file)
            number = 0  # Continue 
            concentration.append(number)

    ####Remove day 4, 5, 11, 12 because of experimental
    concentration = np.delete(concentration, [4,5,11,12], axis=0)

    concentration = pd.Series(concentration, name = "CLL concentration")

    data = pd.concat([viability, concentration], axis=1)

    return data
