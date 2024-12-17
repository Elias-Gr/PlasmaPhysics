
import os
import numpy as np
import json

## all import parameters for "main_optimization"

directory = os.path.dirname(os.path.abspath(__file__))


foldername = "output_2"

USE_ORIGINAL_PARAMETERS = True

##################################################################################
###### ---------------------- Poincare Parameters ----------------------- ########
##################################################################################


NFP = 1 # Number of field periods (be careful if you perturbate a single coils, the number of field periods must be 1)
STELLSYM = False 

PHIPLANES = [0, np.pi/4,  np.pi/2] 

PLOT_PARAMS = {
        "NPOINC": 32,
        "PHIPLANES": PHIPLANES,
        "SAVE": True,
        "SHOW": False,
        "KEEP_RESULT_DATA": False,
        "RLIM": [0.5 * 0.33, 1.4 * 0.33],
        "ZLIM": [-0.42 * 0.33, 0.42 * 0.33]
    }

nLines = 10

RSTART = (np.linspace(0.5, 1.3, nLines)).tolist()
ZSTART = np.zeros(nLines).tolist()

FILEDLINE_PARAMS = {
        "NFP": NFP,
        "STELL_SYM": STELLSYM,
        "R_START": RSTART,
        "Z_START": ZSTART,
        "FOLLOW_TOL": 1e-7,
        "TMAX": 2500,
        "DEGREE": 4,
        "NPLANES": 4

    }

PRESIM = False

###################################################################################
###### ------------------- Optimization Parameters ----------------------- ########
###################################################################################

PRESIM_OPT = False

CUSTOM_CURRENTS = False

OUTER_SURFACE = False

MAXITER = 1234

FIXEDCURRENT = True

r_coil_torus = 0.6

R = 1

r_distance_from_surface = r_coil_torus*0.5


Radius_outer_surface = 3
Radius_outer_surface_z = 2/3 * Radius_outer_surface
Distance_outer_Surface = Radius_outer_surface - r_coil_torus - 1.2
Surface_r = Radius_outer_surface - Distance_outer_Surface
Surface_z = Radius_outer_surface_z - Distance_outer_Surface

nphi = 32
ntheta = 32

CURRENT = 1e5

fourierordercoils = 2

ncoils = 3

WEIGHT_DIST = 1000
WEIGHT_DIST_OUT = 1000
WEIGHT_CURVE = 1
CCDIST_WEIGHT = 1000
BDOTN_WEIGHT = 10000


CURVATURE_THRESHOLD = 8.5
DIST_THRESHOLD = r_distance_from_surface
DIST_THRESHOLD_OUT = Distance_outer_Surface
CCDIST_THRESH = 0.05


######## -------- save parameters ----------- #########

if __name__ == "__main__":

    os.chdir(directory)

    if not os.path.exists(foldername):
        os.mkdir(foldername)

    user_defined_vars = {name: value for name, value in globals().items() if not name.startswith("__")}
    user_defined_vars['os'] = None
    user_defined_vars['np'] = None
    user_defined_vars['json'] = None

    with open(foldername+"/parameters.json", 'w') as file:
        json.dump(user_defined_vars, file, indent = 4)
