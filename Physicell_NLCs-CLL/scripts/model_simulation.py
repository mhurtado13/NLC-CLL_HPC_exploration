import xml.etree.ElementTree as ET
import subprocess
import pandas as pd
import os
import random
from collect_data import collect
from merge_data import merge

def run_model(input_file_path, replicates, *args):                
    
    thread = args[0] #Extract in which thread we are
    values = args[1:]
    errors = []
    tree = ET.parse(input_file_path) #Load xml file
    root = tree.getroot()
    print("Running simulation with parameters " + str(values) + "in node ", str(thread))
    output_folder = "output_" + str(thread)
    param_element = root.find(".//save/folder") #Find the random seed in XML file
    param_element.text = output_folder

    xml_file = "config/configuration_file_" + str(thread) + ".xml"

    os.makedirs(output_folder, exist_ok=True)

    param_behaviors = {'cancer':{'uptake_rate': 0, 'speed': 1, 'transformation_rate': 2},
                    'monocytes':{'speed': 3, 'dead_phagocytosis_rate': 4},
                    'macrophages':{'speed': 5, 'dead_phagocytosis_rate': 6},
                    'NLCs': {'secretion_rate': 7, 'speed': 8, 'dead_phagocytosis_rate': 9},
                    'apoptotic':{'death_rate':10, 'secretion_rate':11}}
        
    for i, celltype in enumerate(param_behaviors.keys()): #i = number of keys name and celltype = cell type
        for param, column in param_behaviors[celltype].items(): #param = parameter name and column = column number
            if celltype == 'cancer' and param == 'uptake_rate':
                param_value = values[column] #Extract each value [i, col_index]
                param_element = root.find(f".//*[@name='{celltype}']//*[@name='anti-apoptotic factor']//{param}") #Find the param name in XML file
                param_element.text = str(param_value)
            elif celltype == 'cancer' and param == 'transformation_rate':
                param_value = values[column] #Extract each value [i, col_index]
                param_element = root.find(f".//*[@name='{celltype}']//{param}/[@name='apoptotic']") #Find the param name in XML file
                param_element.text = str(param_value)
            elif celltype == 'NLCs' and param == 'secretion_rate':
                param_value = values[column] #Extract each value [i, col_index]
                param_element = root.find(f".//*[@name='{celltype}']//*[@name='anti-apoptotic factor']//{param}") #Find the param name in XML file
                param_element.text = str(param_value)
            elif celltype == 'apoptotic' and param == 'death_rate':
                param_value = values[column]
                param_element = root.find(f".//*[@name='{celltype}']//*[@name='apoptosis']//{param}")
                param_element.text = str(param_value)
            elif celltype == 'apoptotic' and param == 'secretion_rate':
                param_value = values[column]
                param_element = root.find(f".//*[@name='{celltype}']//*[@name='debris']//{param}")
                param_element.text = str(param_value)    
            else:
                param_value = values[column] #Extract each value [i, col_index]
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
            with subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE) as proc:
                stdout, stderr = proc.communicate()

            # Check that the Physicell ran successfully
            if proc.returncode != 0:
                print("Error running Physicell")
                print("Physicell error for parameters: " + str(values))
                print(stderr)
                errors.append(values)
                terminate = True
                break #exit loop to avoid running all replicates if there is an error in simulation 

            if terminate == False:
                print("Collecting data in node " + str(thread))
                res = collect(output_folder, xml_file) #We collect the data at each iteration
                data = pd.concat([res, data], axis=1)
        
        loop = False #when the loop finishes exit the while, this is used when terminate == False, meaning that it run succesfully     
          
    if terminate == False:
        viability, concentration = merge(data) #Merge data of replicates 
        print("Physicell simulation for node " + str(thread) + " with parameters " + str(values) + " completed succesfully! :)")
    else:
        viability = pd.Series([0] * 10)
        concentration = pd.Series([0] * 10)
        print("Physicell simulation for node " + str(thread) + " with parameters " + str(values) + " did not run succesfully... completing with 0s")


    return viability, concentration, errors
