
import os
import numpy as np
import json

## all import parameters for "main_optimization"

directory = os.path.dirname(os.path.abspath(__file__))

# Result folder for all your result files and plots
foldername = "test_Radius_0.55" 

# Vmec file of plasma surface to use for optmization
plasma_surface_filename = "input" 

USE_ORIGINAL_PARAMETERS = False

##################################################################################
###### ---------------------- Poincare Parameters ----------------------- ########
##################################################################################

# Number of field periods (be careful if you perturbate a single coils, the number of field periods must be 1)

NFP = 1 

STELLSYM = False 

PHIPLANES = [np.pi/8, np.pi/4,  np.pi/2] 

PLOT_PARAMS = {
        "NPOINC": 32,
        "PHIPLANES": PHIPLANES,
        "SAVE": True,
        "SHOW": False,
        "KEEP_RESULT_DATA": False,
        "RLIM": [1, 1.3],
        "ZLIM": [0,0]#[-0.42 * 0.33, 0.42 * 0.33]
    }

# define starting points of the space curves for the poincare plots

nLines = 30

RSTART = (np.linspace(0.5, 1.3, nLines)).tolist()
ZSTART = np.zeros(nLines).tolist()

FILEDLINE_PARAMS = {
        "NFP": NFP,
        "STELL_SYM": STELLSYM,
        "R_START": RSTART,
        "Z_START": ZSTART,
        "FOLLOW_TOL": 1e-10,
        "TMAX": 25000,
        "DEGREE": 4,
        "NPLANES": 4

    }

# use the pre-optimized from vtu file


PRESIM = False

###################################################################################
###### ------------------- Optimization Parameters ----------------------- ########
###################################################################################


PRESIM_OPT = False

# maximum iterations of minimization algorithm:

MAXITER = 1234

# minor radius of torus to fix the coils onto:

r_coil_torus = 0.6

# major radius of torus:

R = 1

# number of angle points for plasma surface

nphi = 32
ntheta = 32

#  starting current of the coils:

CURRENT = 1e5

FIXEDCURRENT = True

CUSTOM_CURRENTS = False

# fourierorder of the coils used for optimization (the more the more complex)

fourierordercoils = 2

# number of coils per quarter:

ncoils = 3

### weights for the different cost function objects:

DISTANCE_TORUS = 1000
DISTANCE_OUTER_SURFACE = 1000
CURVATURE = 1
CURVECURVEDISTANCE = 1000
FLUX = 10000


# threshholds for different cost function measures

CURVATURE_THRESH = 8.5

CURVECURVEDIST_THRESH = 0.05


######## -------- save parameters ----------- #########

if __name__ == "__main__":

    os.chdir(directory)

    if not os.path.exists(foldername):
        os.mkdir(foldername)

    user_defined_vars = {name: value for name, value in globals().items() if not (name.startswith("__") or "np" or "os" or "json")}


    with open(foldername+"/parameters.json", 'w') as file:
        json.dump(user_defined_vars, file, indent = 4)
