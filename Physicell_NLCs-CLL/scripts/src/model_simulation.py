import xml.etree.ElementTree as ET
import subprocess
import pandas as pd
import os
import shutil
import random
from collect_data import collect
from merge_data import merge

def simulate_model(input_file_path, replicates, node, process, *args):                
    
    values = args
    errors = []
    tree = ET.parse(input_file_path) #Load original xml file
    root = tree.getroot()
    print("Running simulation with parameters " + str(values) + " in task number " + str(process) + " using node " + str(node) + "\n")

    #Creating folder for each node
    output_node = "output/output_node_" + str(node)
    os.makedirs(output_node, exist_ok=True)
    
    #Creating subfolder for each process inside each node
    output_folder = output_node + "/output_" + str(process)
    os.makedirs(output_folder, exist_ok=True)

    #Parsing output folder in xml file 
    param_element = root.find(".//save/folder") #Find the random seed in XML file
    param_element.text = output_folder

    #Creating subfolder for each node inside config
    config_node = "config/config_node_" + str(node)
    os.makedirs(config_node, exist_ok=True)

    xml_file = config_node + "/configuration_file_" + str(process) + ".xml"

    ### All parameters
    ''' param_behaviors = {
        'cancer':{'cell_cell_repulsion_strength': 0, 
                  'speed': 1, 
                  'secretion_rate': {'column': 2, 'substrate': 'cytokines'}, 
                  'uptake_rate': {'column': 3, 'substrate': 'antiapoptotic'},
                  'transformation_rate': {'column':4, 'cell_type': 'apoptotic'}},

        'monocytes':{'speed': 5, 
                     'uptake_rate': {'column': 6, 'substrate': 'cytokines'},
                     'uptake_rate': {'column': 7, 'substrate': 'stress'},
                     'phagocytosis_rate': {'column':8, 'cell_type': 'apoptotic'},
                     'phagocytosis_rate': {'column':9, 'cell_type': 'dead'},
                     'transformation_rate': {'column':10, 'cell_type': 'macrophages'},
                     'transformation_rate': {'column':11, 'cell_type': 'NLCs'}},

        'macrophages':{'speed': 12, 
                       'uptake_rate': {'column': 13, 'substrate': 'cytokines'},
                       'uptake_rate': {'column': 14, 'substrate': 'stress'},
                       'attack_rate': {'column': 15, 'substrate': 'cancer'},
                       'damage_rate': 16, 
                       'phagocytosis_rate': {'column': 17, 'cell_type': 'apoptotic'},
                       'phagocytosis_rate': {'column': 18, 'cell_type': 'dead'},
                       'transformation_rate': {'column': 19, 'cell_type': 'NLCs'}},

        'NLCs': {'cell_cell_adhesion_strength': 20, 
                 'attachment_rate': 21, 
                 'detachment_rate': 22,
                 'speed': 23, 
                 'secretion_rate': {'column': 24, 'substrate': 'antiapoptotic'},
                 'uptake_rate': {'column': 25, 'substrate': 'stress'},
                 'phagocytosis_rate': {'column': 26, 'cell_type': 'apoptotic'},
                 'phagocytosis_rate': {'column': 27, 'cell_type': 'dead'}},

        'apoptotic':{'speed': 28, 
                     'secretion_rate': {'column': 29, 'substrate': 'cytokines'},
                     'secretion_rate': {'column': 30, 'substrate': 'stress'},
                     'transformation_rate': {'column': 31, 'cell_type': 'dead'}},

        'dead':{'secretion_rate': {'column': 32, 'substrate': 'stress'}}
    } '''

    ### Parameters selected after Exploration

    '''param_behaviors = {
        'cancer':{'cell_cell_repulsion_strength': 0, 
                  'speed': 1},

        'monocytes':{'phagocytosis_rate': {'column':2, 'cell_type': 'dead'},
                     'transformation_rate': {'column':3, 'cell_type': 'NLCs'}},

        'macrophages':{'uptake_rate': {'column': 4, 'substrate': 'cytokines'},
                       'attack_rate': {'column': 5, 'substrate': 'cancer'},
                       'phagocytosis_rate': {'column': 6, 'cell_type': 'apoptotic'},
                       'phagocytosis_rate': {'column': 7, 'cell_type': 'dead'},
                       'transformation_rate': {'column': 8, 'cell_type': 'NLCs'}},

        'NLCs': {'cell_cell_adhesion_strength': 9, 
                 'attachment_rate': 10, 
                 'phagocytosis_rate': {'column': 11, 'cell_type': 'dead'}},

        'apoptotic':{'speed': 12, 
                     'transformation_rate': {'column': 13, 'cell_type': 'dead'}}
    }'''

    ### Parameters selected after SA

    param_behaviors = {
        'cancer':{'cell_cell_repulsion_strength': 0, 
                  'speed': 1},

        'monocytes':{'phagocytosis_rate': {'column':2, 'cell_type': 'dead'}},

        'macrophages':{'uptake_rate': {'column': 3, 'substrate': 'cytokines'},
                       'attack_rate': {'column': 4, 'substrate': 'cancer'}},

        'NLCs': {'cell_cell_adhesion_strength': 5, 
                 'phagocytosis_rate': {'column': 6, 'cell_type': 'dead'}},

        'apoptotic':{'speed': 7, 
                     'transformation_rate': {'column': 8, 'cell_type': 'dead'}}
    }
                        
    for _, celltype in enumerate(param_behaviors.keys()): #i = number of keys name and celltype = cell type
        for param, column_info in param_behaviors[celltype].items(): #param = parameter name and column = column number
            if isinstance(column_info, dict):
                column = column_info['column']
                extra = column_info.get('cell_type') or column_info.get('substrate')
            else:
                column = column_info
                extra = None

            param_value = values[column] #Extract each value [i, col_index]

            if param in ['uptake_rate', 'secretion_rate'] and extra:
                param_element = root.find(f".//*[@name='{celltype}']//*[@name='{extra}']//{param}") #Find the param name in XML file
            elif extra:
                param_element = root.find(f".//*[@name='{celltype}']//{param}/[@name='{extra}']") #Find the param name in XML file
            else:
                param_element = root.find(f".//*[@name='{celltype}']//{param}") #Find the param name in XML file

            param_element.text = str(param_value)

    # Define the command to call your C++ software with the updated XML as input
    command = ["./project", xml_file]
    data = pd.DataFrame()        
    terminate = False
    loop = True
    while not terminate and loop:
        for i in range(replicates): #replicates is for bootstrapping, we run the simulation with updated value # (replicates) times
            # Random seed for each simulation
            param_element = root.find(".//random_seed") #Find the random seed in XML file
            param_element.text = str(random.randint(0,4294967295))

            # Write the updated XML to a string
            updated_xml_str = ET.tostring(root, encoding="unicode", method="xml")

            with open(xml_file, "w") as file:
                file.write(updated_xml_str)

            # Call the C++ software using subprocess
            print("Running model in task " + str(process) + "\n")
            with subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE) as proc:
                stdout, stderr = proc.communicate()

            # Check that the Physicell ran successfully
            if proc.returncode != 0:
                print("Error running Physicell\n")
                print(stderr)
                errors.append(values)
                terminate = True
                break #exit loop to avoid running all replicates if there is an error in simulation 

            if terminate == False:
                print("Collecting data in task " + str(process) + "\n")
                res = collect(output_folder, xml_file, node) #We collect the data at each iteration
                data = pd.concat([res, data], axis=1)
        
        loop = False #when the loop finishes exit the while, this is used when terminate == False, meaning that it run succesfully     
          
    if terminate == False:
        viability, concentration = merge(data) #Merge data of replicates 
        print("Physicell simulation for task " + str(process) + " with parameters " + str(values) + " in node " + str(node) + " completed succesfully! :)\n")
    else:
        viability = pd.Series([0] * 10)
        concentration = pd.Series([0] * 10)
        print("Physicell simulation for task " + str(process) + " with parameters " + str(values) + " in node " + str(node) + " did not run succesfully... completing with 0s\n")

    #Clean memory
    shutil.rmtree(output_folder) #Remove output folder for each task when it is done 
    os.remove(xml_file) #Remove xml file for each task when it is done

    return viability, concentration, errors
