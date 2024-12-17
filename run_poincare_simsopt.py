import simsopt.field
from simsopt.geo import SurfaceRZFourier, plot
import runPoincare as rp
import pertubateCoils as pc
import numpy as np
import simsopt
import os
import json
from simsopt._core import load
from simsopt.field import BiotSavart, coils_via_symmetries
from parameters import *


####### DESIGN INPUTS ######

os.chdir(directory)

NFP = 1 # Number of field periods (be careful if you perturbate a single coils, the number of field periods must be 1)
STELLSYM = False 
LOG_FIELDLINES = True       # If True, print the output of xfiledlines to the console

input_file_path = foldername + '/bs.json'
surface_filename = 'input'

output_dir = os.path.join(foldername) 
if not os.path.exists(output_dir):
    os.makedirs(output_dir)


####### BASELINE ######
# Set to True to run without perturbations
RUN_BASELINE = True

####### PERTURBATIONS ######
# Set to True to save BioSavar opject of perturbaated coils
SAVE_PERTURBATIONS = False

# Set to True to run with shifted coils
RUN_SHIFTED_COILS = False
# Set array of coils to shift
SHIFTED_COILS = [0]
# Set the shifts to apply
SHIFTS = np.arange(0.03, 0.051, 0.01) #np.arange(0.01, 0.1, 0.01) * 0.33#np.linspace(0.0001, 0.01, 10)

# Set to True to run with perturbated current
RUN_PERTURBATED_CURRENT = False
# Set array of coils to perturbate
PERTURBATED_COILS = [0,3,6,9]
# Set the perturbation amplitude in percentage (0.1 = 10%)
PERTURBATION_AMPLITUDE = np.linspace(0.0001, 0.01, 5)


#---------------------------------#
# load the file acording to your input file. Outcommand the one you are not using

# load the input file option 1
# surfaces, B_axis, coils = simsopt.load(input_file_path)
# bs = simsopt.field.BiotSavart(coils)

# load the input file option 2
#

surface = SurfaceRZFourier.from_vmec_input(surface_filename, range="full torus", nphi=64, ntheta=64)

if PRESIM:
    coils_presim = load(f'serial0021326.json')[1]
    base_curves = [x.curve for x in coils_presim]
    base_currents = [x.current for x in coils_presim]
    coils = coils_via_symmetries(base_curves, base_currents, surface.nfp, True)
    
    bs = BiotSavart(coils)
else:
    bs = simsopt.load(input_file_path)



coils = bs.coils
for i, c in enumerate(coils):
    print(f'Coil {i} Current: {c.current.get_value()}')




FILEDLINE_PARAMS["EXTCUR"] = [c.current.get_value() for c in coils]

#---------------------------------#
if RUN_BASELINE:
    print("Running baseline...")
    # Check if the directory exists and create it if it doesn't
    output_dir = os.path.join(foldername, 'BASELINE')
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    variant = 'BASELINE'
    rp.run_poincare_simsopt(coils, variant, output_dir, surface, SAVE_PERTURBATIONS, FILEDLINE_PARAMS)
 

if RUN_SHIFTED_COILS:
    print("Running shifted coils...")
    # Check if the directory exists and create it if it doesn't
    output_dir = os.path.join(foldername, 'SHIFTS')
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for shift in SHIFTS:
        variant = f'{shift:.4f}_Shift'
        shifted_coils = coils.copy()
        for i in SHIFTED_COILS:
            trans_vec = np.array([shift, 0, 0])
            shifted_coil = pc.translate_coil(coils[i], trans_vec)
            shifted_coils[i] = shifted_coil
        
        # plot([surface] + shifted_coils, engine="mayavi", close=True)
        rp.run_poincare_simsopt(shifted_coils, variant, output_dir, surface, SAVE_PERTURBATIONS, FILEDLINE_PARAMS)

if RUN_PERTURBATED_CURRENT:
    print("Running perturbated current...")
    # Check if the directory exists and create it if it doesn't
    output_dir = os.path.join(foldername, 'SIMSOPT', 'CURRENT_PERTURBATION')
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for p in PERTURBATION_AMPLITUDE:
        variant = f'{p:.3f}_Current_Perturbation'
        pertupated_coils = coils.copy()   
        for i in PERTURBATED_COILS:
            current = coils[i].current
            print(f'unpertubated current: {current.get_value()}')
            p_current  = current * (1+p)
            pertupated_coils[i] = simsopt.field.Coil(coils[i].curve, p_current)
            print(f'pertubated current: {pertupated_coils[i].current.get_value()}')
        
        FILEDLINE_PARAMS["EXTCUR"] = [c.current.get_value() for c in pertupated_coils]
        rp.run_poincare_simsopt(pertupated_coils, variant, output_dir, surface, SAVE_PERTURBATIONS, FILEDLINE_PARAMS)

