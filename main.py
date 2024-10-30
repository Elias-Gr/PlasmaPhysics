""" #%%
from simsopt.geo import CurveXYZFourier

curve = CurveXYZFourier(300,1)
curve.set_dofs([0,0,0,0,1,0,0,0,1])
curve.plot() """

import matplotlib.pyplot as plt
import numpy as np
from simsopt._core import load
from simsopt.geo import plot, MeanSquaredCurvature
from simsopt.geo import SurfaceRZFourier, create_equally_spaced_curves, \
    CurveLength, curves_to_vtk
from simsopt.field import Current, coils_via_symmetries, BiotSavart
from simsopt.objectives import SquaredFlux, QuadraticPenalty
from scipy.optimize import minimize

# replace "NAME_OF_FILE_YOU_DOWNLOADED" with the name you gave the file
res = load(f'serial0021326.json')

ncoils = 4

MAXITER = 300

nphi = 64
ntheta = 64
s = SurfaceRZFourier.from_vmec_input('input', range="full torus", nphi=nphi, ntheta=ntheta)

coils = res[1]
surface = res[0]

print(np.size(coils))

##creating coils and currents:

base_curves = create_equally_spaced_curves(ncoils, s.nfp, stellsym=True, R0=1, R1=0.6, order=5)
base_currents = [Current(1.0) * 1e5 for i in range(ncoils)]
#coils = coils_via_symmetries(base_curves, base_currents, s.nfp, True)

##calculating magnetic field from coils:

bs = BiotSavart(coils)
bs.set_points(s.gamma().reshape((-1, 3)))

##setting up parameters for penalty function:

Jf = SquaredFlux(s, bs)
Jls = [CurveLength(c) for c in base_curves]

LENGTH_WEIGHT=0
LENGTH_TARGET = 18.0
C_WEIGHT = 1
MSC_THRESHOLD = 10 #1e-6

Jl = QuadraticPenalty(sum(Jls), LENGTH_TARGET, "max")
Jmscs = [MeanSquaredCurvature(c) for c in base_curves]
Jc = sum(QuadraticPenalty(J, MSC_THRESHOLD) for J in Jmscs)

## penalty function: (change penalties fitting to your problem)
## (here it's length of coils curves. use radius for 1st problem i guess)

JF = Jf + LENGTH_WEIGHT * Jl + C_WEIGHT * Jc


B_dot_n = np.sum(bs.B().reshape((nphi, ntheta, 3)) * s.unitnormal(), axis=2)

print('Initial max|B dot n|:', np.max(np.abs(B_dot_n)))


plot(coils + [s], engine="plotly", close=True)

def fun(dofs):
    JF.x = dofs
    return JF.J(), JF.dJ()

print("""
################################################################################
### Perform a Taylor test ######################################################
################################################################################
""")
f = fun
dofs = JF.x
np.random.seed(1)
h = np.random.uniform(size=dofs.shape)
J0, dJ0 = f(dofs)
dJh = sum(dJ0 * h)
for eps in [1e-3, 1e-4, 1e-5, 1e-6, 1e-7]:
    J1, _ = f(dofs + eps*h)
    J2, _ = f(dofs - eps*h)
    print("err", (J1-J2)/(2*eps) - dJh)

print("""
################################################################################
### Run the optimisation #######################################################
################################################################################
""")

res = minimize(fun, dofs, jac=True, method='L-BFGS-B',
               options={'maxiter': MAXITER, 'maxcor': 500, 'iprint': 5}, tol=1e-15)

B_dot_n = np.sum(bs.B().reshape((nphi, ntheta, 3)) * s.unitnormal(), axis=2)
print('Final max|B dot n|:', np.max(np.abs(B_dot_n)))
total_curve_length = sum(func.J() for func in Jls)
print("Sum of lengths of base curves:", total_curve_length)


plot(coils + [s], engine="plotly", close=True)
### minimization implementation