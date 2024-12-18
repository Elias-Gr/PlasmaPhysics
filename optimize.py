import sys
import os
import subprocess
import json

directory = os.path.dirname(os.path.abspath(__file__))

os.chdir(directory)


def change_parameter(parameter, new_value, original=False):
    x = False
    if original:
        # File to modify
        file_path = sys.argv[1]+"/parameters.json"
        # Check if the file exists
        if os.path.exists(file_path):
            # Read the file and replace the value
            with open(file_path,'r') as data:
                params_dict = json.load(data)
            typ = type(params_dict[parameter])
            params_dict[parameter] = typ(new_value)  
            with open(file_path,'w') as file:
                json.dump(params_dict, file, indent=4)
            print(f'dict entry {parameter} in {file_path} changed to {new_value}')
            

    else:
        # File to modify
        file_path = "config.py"
        # Check if the file exists
        if os.path.exists(file_path):
            # Read the file and replace the value
            with open(file_path, "r") as file:
                lines = file.readlines()

            with open(file_path, "w") as file:
                for line in lines:
                    if line.startswith(parameter+" ="):
                        file.write(parameter+f" = {new_value}\n")  # Replace the value
                        print(f"Successfully updated {parameter} in {file_path} to '{new_value}'.")
                        x = True
                    else:
                        file.write(line)  # Keep other lines unchanged
        
        else:
            print(f"Error: {file_path} does not exist.")
        if not x:
            print(f"Error: {parameter} is not part of config.py parameters")


try:

    output_folder = sys.argv[1]
except:
    subprocess.run(['python', 'config.py'])
    subprocess.run(['python', 'main_optimization.py'])
    subprocess.run(['python', 'run_poincare_simsopt.py'])
    exit()


change_parameter('foldername', f'\"{output_folder}\"')

if '-orp' in sys.argv:
    change_parameter('USE_ORIGINAL_PARAMETERS',True)
    original_parameters = True
else:
    change_parameter('USE_ORIGINAL_PARAMETERS',False)
    original_parameters = False

if '-c' in sys.argv:
    for i, arg in enumerate(sys.argv):
        if arg == '-c':
            parameter = sys.argv[i+1]
            value = sys.argv[i+2]
            change_parameter(parameter,value,original=original_parameters)


if '-o' in sys.argv: 
    if not original_parameters:
        subprocess.run(['python', 'config.py'])
    subprocess.run(['python', 'main_optimization.py'])
elif '-p' in sys.argv:
    subprocess.run(['python', 'run_poincare_simsopt.py'])
elif '-n' in sys.argv:
    pass
else:
    if not original_parameters:
        subprocess.run(['python', 'config.py'])
    subprocess.run(['python', 'main_optimization.py'])
    subprocess.run(['python', 'run_poincare_simsopt.py'])


