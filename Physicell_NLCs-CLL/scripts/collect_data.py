import numpy as np
import pandas as pd
import pcdl

#config_file = 'config/NLC_CLL.xml'
def collect(dir_output, config_file):
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

    positions = []
    for days in range(0,14):
        hours = 24*days
        positions.append(hours)

    alive = [CLL_initial]
    dead = [0]
    apoptotic = [apoptotic_initial]
    for i in range(1, len(positions)):
        step = timesteps[positions[i]].get_cell_df(states=1)
        number_alive = len(step[(step['cell_type']=='cancer')&(step['dead']==False)])
        number_apoptotic = len(step[(step['cell_type']=='apoptotic')&(step['dead']==False)])
        number_dead = len(step[step['dead']==True])
        alive.append(number_alive)
        dead.append(number_dead)
        apoptotic.append(number_apoptotic)

    CLL_alive = pd.Series(alive, name="Cells_alive")
    CLL_apoptotic = pd.Series(apoptotic, name = "Cells_apoptotic")
    CLL_dead = pd.Series(dead, name = "Cells_dead")

    #viability at time t =  CLL alive at time t / (CLL alive + CLL apoptotic + CLL dead) at time t
    viability = []
    for i in range(len(CLL_alive)):
        number = (CLL_alive[i]/(CLL_alive[i]+CLL_apoptotic[i]+CLL_dead[i]))*100
        viability.append(number)

    ####Remove day 4, 5, 11, 12 because of experimental
    viability = np.delete(viability, [4,5,11,12], axis=0)

    viability = pd.Series(viability, name = "CLL viability")

    #concentration at time t =  CLL alive at time t / (CLL initial)
    concentration = []
    for i in range(len(CLL_alive)):
        number = (CLL_alive[i]/CLL_initial)*100
        concentration.append(number)

    ####Remove day 4, 5, 11, 12 because of experimental
    concentration = np.delete(concentration, [4,5,11,12], axis=0)

    concentration = pd.Series(concentration, name = "CLL concentration")

    data = pd.concat([viability, concentration], axis=1)

    return data
    #file_csv = 'data_output/data.csv'

    #If the file already exists
    #if os.path.exists(file_csv):
    #    old_data = pd.read_csv(file_csv)
    #    new_data = pd.concat([old_data, df], axis=1)
    #    new_data.to_csv(file_csv, index=False, header=True)
    #else:
    #    df.to_csv(file_csv, index=False, header=True)