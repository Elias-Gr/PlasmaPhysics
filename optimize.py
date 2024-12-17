import sys
import os
import subprocess

directory = os.path.dirname(os.path.abspath(__file__))

os.chdir(directory)


output_folder = sys[1]


# File to modify
file_path = "config.py"

def change_parameter(parameter, new_value):
# Check if the file exists
    if os.path.exists(file_path):
        # Read the file and replace the value
        with open(file_path, "r") as file:
            lines = file.readlines()

        with open(file_path, "w") as file:
            for line in lines:
                if line.startswith(parameter+" ="):
                    file.write(parameter+f" = \"{new_value}\"\n")  # Replace the value
                else:
                    file.write(line)  # Keep other lines unchanged

        print(f"Successfully updated {parameter} in {file_path} to '{new_value}'.")
    else:
        print(f"Error: {file_path} does not exist.")

change_parameter('foldername',output_folder)

exit


subprocess.run(['python', 'main_optimization.py'])
subprocess.run(['python', 'run_poincare_simsopt.py'])
