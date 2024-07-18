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
    print("Running simulation with parameters " + str(values) + " in task number " + str(process) + " using node " + str(node))
    
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

    param_behaviors = {'cancer':{'uptake_rate': 0, 'secretion_rate': 1, 'speed': 2, 'transformation_rate': 3, 'relative_maximum_adhesion_distance': 4, 'cell_adhesion_affinity': 5},
                    'monocytes':{'speed': 6, 'transformation_rate': 7, 'phagocytosis_rate': 8}, 
                    'macrophages':{'speed': 9, 'phagocytosis_rate': 10, 'attack_rate': 11, 'relative_maximum_adhesion_distance': 12, 'cell_adhesion_affinity': 13},
                    'NLCs': {'speed': 14, 'phagocytosis_rate': 15},
                    'apoptotic':{'secretion_rate':16, 'speed':17, 'transformation_rate_apoptotic': 18},
                    'dead':{'secretion_rate_dead': 19}}
        
    for i, celltype in enumerate(param_behaviors.keys()): #i = number of keys name and celltype = cell type
        for param, column in param_behaviors[celltype].items(): #param = parameter name and column = column number
            if celltype == 'cancer' and param == 'uptake_rate':
                param_value = values[column] #Extract each value [i, col_index]
                param_element = root.find(f".//*[@name='{celltype}']//*[@name='anti-apoptotic factor']//{param}") #Find the param name in XML file
                param_element.text = str(param_value)
            if celltype == 'cancer' and param == 'secretion_rate':
                param_value = values[column] #Extract each value [i, col_index]
                param_element = root.find(f".//*[@name='{celltype}']//*[@name='cancer-signal']//{param}") #Find the param name in XML file
                param_element.text = str(param_value)
            elif celltype == 'cancer' and param == 'transformation_rate':
                param_value = values[column] #Extract each value [i, col_index]
                param_element = root.find(f".//*[@name='{celltype}']//{param}/[@name='apoptotic']") #Find the param name in XML file
                param_element.text = str(param_value)
            elif celltype == 'cancer' and param == 'cell_adhesion_affinity':
                param_value = values[column] #Extract each value [i, col_index]
                param_element = root.find(f".//*[@name='{celltype}']//{param}/[@name='NLCs']") #Find the param name in XML file
                param_element.text = str(param_value)
            elif celltype == 'monocytes' and param == 'transformation_rate':
                param_value = values[column] #Extract each value [i, col_index]
                param_element = root.find(f".//*[@name='{celltype}']//{param}/[@name='macrophages']") #Find the param name in XML file
                param_element.text = str(param_value)
            elif celltype == 'monocytes' and param == 'phagocytosis_rate':
                param_value = values[column] #Extract each value [i, col_index]
                param_element = root.find(f".//*[@name='{celltype}']//{param}/[@name='apoptotic']") #Find the param name in XML file
                param_element.text = str(param_value)
            elif celltype == 'macrophages' and param == 'phagocytosis_rate':
                param_value = values[column] #Extract each value [i, col_index]
                param_element = root.find(f".//*[@name='{celltype}']//{param}/[@name='apoptotic']") #Find the param name in XML file
                param_element.text = str(param_value)
            elif celltype == 'macrophages' and param == 'attack_rate':
                param_value = values[column] #Extract each value [i, col_index]
                param_element = root.find(f".//*[@name='{celltype}']//{param}/[@name='cancer']") #Find the param name in XML file
                param_element.text = str(param_value)
            elif celltype == 'macrophages' and param == 'cell_adhesion_affinity':
                param_value = values[column] #Extract each value [i, col_index]
                param_element = root.find(f".//*[@name='{celltype}']//{param}/[@name='cancer']") #Find the param name in XML file
                param_element.text = str(param_value)                    
            elif celltype == 'NLCs' and param == 'phagocytosis_rate':
                param_value = values[column] #Extract each value [i, col_index]
                param_element = root.find(f".//*[@name='{celltype}']//*[@name='apoptotic']//{param}") #Find the param name in XML file
                param_element.text = str(param_value)
            elif celltype == 'apoptotic' and param == 'secretion_rate':
                param_value = values[column]
                param_element = root.find(f".//*[@name='{celltype}']//*[@name='dead substrate']//{param}")
                param_element.text = str(param_value)
            elif celltype == 'apoptotic' and param == 'transformation_rate':
                param_value = values[column] #Extract each value [i, col_index]
                param_element = root.find(f".//*[@name='{celltype}']//{param}/[@name='dead']") #Find the param name in XML file
                param_element.text = str(param_value)
            elif celltype == 'dead' and param == 'secretion_rate':
                param_value = values[column]
                param_element = root.find(f".//*[@name='{celltype}']//*[@name='dead substrate']//{param}")
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
                print(stderr)
                errors.append(values)
                terminate = True
                break #exit loop to avoid running all replicates if there is an error in simulation 

            if terminate == False:
                print("Collecting data in task " + str(process))
                res = collect(output_folder, xml_file) #We collect the data at each iteration
                data = pd.concat([res, data], axis=1)
        
        loop = False #when the loop finishes exit the while, this is used when terminate == False, meaning that it run succesfully     
          
    if terminate == False:
        viability, concentration = merge(data) #Merge data of replicates 
        print("Physicell simulation for task " + str(process) + " with parameters " + str(values) + " in node " + str(node) + " completed succesfully! :)")
    else:
        viability = pd.Series([0] * 10)
        concentration = pd.Series([0] * 10)
        print("Physicell simulation for task " + str(process) + " with parameters " + str(values) + " in node " + str(node) + " did not run succesfully... completing with 0s")

    return viability, concentration, errors
