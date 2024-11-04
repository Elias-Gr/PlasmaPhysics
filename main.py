""" #%%
from simsopt.geo import CurveXYZFourier

curve = CurveXYZFourier(300,1)
curve.set_dofs([0,0,0,0,1,0,0,0,1])
curve.plot() """

import matplotlib.pyplot as plt
import numpy as np
from simsopt._core import load
from simsopt.geo import plot, MeanSquaredCurvature, CurveSurfaceDistance, CurveCurveDistance, LpCurveCurvature
from simsopt.geo import SurfaceRZFourier, create_equally_spaced_curves, \
    CurveLength, curves_to_vtk
from simsopt.field import Current, coils_via_symmetries, BiotSavart
from simsopt.objectives import SquaredFlux, QuadraticPenalty
from scipy.optimize import minimize

res = load(f'serial0021326.json')

MAXITER = 600

USE_PRESIM = False




nphi = 48
ntheta = 48
s = SurfaceRZFourier.from_vmec_input('input', range="full torus", nphi=nphi, ntheta=ntheta)

##creating coils and currents:
if USE_PRESIM:
    coils_presim = res[1]
    surface = res[0][-1]
    base_curves = [x.curve for x in coils_presim]
    base_currents = [x.current for x in coils_presim]
    base_currents[0].fix_all()
else:
    ncoils = 4
    base_curves = create_equally_spaced_curves(ncoils, s.nfp, stellsym=True, R0=1, R1=0.6, order=5)
    base_currents = [Current(1.0) * 1e5 for i in range(ncoils)]
    base_currents[0].fix_all()
   
   


coils = coils_via_symmetries(base_curves, base_currents, s.nfp, True)

plot(coils + [s], engine="plotly", close=True)



##calculating magnetic field from coils:

bs = BiotSavart(coils)
bs.set_points(s.gamma().reshape((-1, 3)))

##setting up parameters for penalty function:

Jf = SquaredFlux(s, bs)
Jls = [CurveLength(c) for c in base_curves]

LENGTH_WEIGHT = 0.01
C_WEIGHT = 0.05
DIST_WEIGHT = 1
CCD_WEIGHT = 1

LENGTH_TARGET = 13.5
MSC_THRESHOLD = 10
CURVATURE_THRESHOLD = 10
DIST_THRESHOLD = 0.1
CCDIST_TH = 0.1
THRESH_CURV = 1

Jl = QuadraticPenalty(sum(Jls), LENGTH_TARGET, "max")
Jmscs = [MeanSquaredCurvature(c) for c in base_curves]
Jcc1 = [LpCurveCurvature(c, 1) for c in base_curves]
Jcc = sum(QuadraticPenalty(J, CURVATURE_THRESHOLD, 'max') for J in Jcc1)
Jc = sum(QuadraticPenalty(J, MSC_THRESHOLD, 'max') for J in Jmscs)


Jd = CurveSurfaceDistance(base_curves, s, DIST_THRESHOLD)

Jccd = CurveCurveDistance(base_curves, CCDIST_TH)

## penalty function: (change penalties fitting to your problem)
## (here it's length of coils curves. use radius for 1st problem i guess)

JF = Jf + Jcc + CCD_WEIGHT * Jccd#+ LENGTH_WEIGHT * Jl + C_WEIGHT * Jc + DIST_WEIGHT * Jd 


B_dot_n = np.sum(bs.B().reshape((nphi, ntheta, 3)) * s.unitnormal(), axis=2)

print('Initial max|B dot n|:', np.max(np.abs(B_dot_n)))




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

B_dot_n = np.sum(bs.B().reshape((nphi, ntheta, 3)) * s.unitnormal(), axis=2)
B_dot_n_max = np.max(np.abs(B_dot_n)).copy()

total_curve_length_start = sum(func.J() for func in Jls).copy()


res = minimize(fun, dofs, jac=True, method='L-BFGS-B',
               options={'maxiter': MAXITER, 'maxcor': 500, 'iprint': 5}, tol=1e-15)

B_dot_n = np.sum(bs.B().reshape((nphi, ntheta, 3)) * s.unitnormal(), axis=2)
print('starting max|B dot n|:', B_dot_n_max)
print('Final max|B dot n|:', np.max(np.abs(B_dot_n)))
total_curve_length = sum(func.J() for func in Jls)
print("starting Sum of lengths of base curves:", total_curve_length_start)
print("Sum of lengths of base curves:", total_curve_length)


plot(coils + [s], engine="plotly", close=True)

