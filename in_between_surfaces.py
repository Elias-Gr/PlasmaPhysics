import matplotlib.pyplot as plt
import numpy as np
from simsopt._core import load
from simsopt.geo import plot, MeanSquaredCurvature, CurveSurfaceDistance, CurveCurveDistance, LpCurveCurvature
from simsopt.geo import SurfaceRZFourier, create_equally_spaced_curves, \
    CurveLength, curves_to_vtk
from simsopt.field import Current, coils_via_symmetries, BiotSavart
from simsopt.objectives import SquaredFlux, QuadraticPenalty
from scipy.optimize import minimize
import simsoptpp as sopp
from scipy.optimize import OptimizeResult

surface_coils = SurfaceRZFourier()
surface_two = SurfaceRZFourier()
surface_plot_coils = SurfaceRZFourier()

MAXITER = 500

Radius = 0.6 #parameter search

distance_from_surfaces = 0.3

surface_coils.set_rc(0,0,1)
surface_coils.set_rc(1,0,Radius-distance_from_surfaces)
surface_coils.set_zs(1,0,Radius-distance_from_surfaces)

surface_two.set_rc(0,0,1)
surface_two.set_rc(1,0,Radius+distance_from_surfaces)
surface_two.set_zs(1,0,Radius+distance_from_surfaces)

surface_plot_coils.set_rc(0,0,1)
surface_plot_coils.set_rc(1,0,Radius)
surface_plot_coils.set_zs(1,0,Radius)

plot([surface_coils], engine='plotly')
plot([surface_two], engine = 'plotly')

nphi = 64
ntheta = 64
s = SurfaceRZFourier.from_vmec_input('input', range="full torus", nphi=nphi, ntheta=ntheta)




ncoils = 3
base_curves = create_equally_spaced_curves(ncoils, s.nfp, stellsym=True, R0=1, R1=Radius, order=3)
base_currents = [Current(1.0) * 1e2 for i in range(ncoils)]
#base_currents[0].fix_all()

coils = coils_via_symmetries(base_curves, base_currents, s.nfp, True)

#plot([surface_coils] + coils +[s], engine="plotly", close=True)

bs = BiotSavart(coils)
bs.set_points(s.gamma().reshape((-1, 3)))

##setting up parameters for penalty function:

CURVATURE_THRESHOLD = 8.5
DIST_THRESHOLD = distance_from_surfaces
CCDIST_THRESH = 0.05

Jf = SquaredFlux(s, bs, definition = 'quadratic flux')

Jcc1 = [LpCurveCurvature(c, 1) for c in base_curves]
Jcc = sum(QuadraticPenalty(J, CURVATURE_THRESHOLD, 'max') for J in Jcc1)


#import pdb; pdb.set_trace()

Jdist_1 = CurveSurfaceDistance(base_curves, surface_coils, DIST_THRESHOLD)
Jdist_2 = CurveSurfaceDistance(base_curves, surface_two, DIST_THRESHOLD)
Jccdist = CurveCurveDistance(base_curves, CCDIST_THRESH)

WEIGHT_DIST = 110
WEIGHT_CURVE = 3
CCDIST_WEIGHT = 100

JF = 1000*Jf + WEIGHT_CURVE * Jcc + WEIGHT_DIST * Jdist_1 + WEIGHT_DIST * Jdist_2 + CCDIST_WEIGHT * Jccdist


##### minimization #####


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

#total_curve_length_start = sum(func.J() for func in Jls).copy()

call_result = []

def callback(intermediate_result: OptimizeResult):
    call_result.append(intermediate_result['fun'])

res = minimize(fun, dofs, jac=True, method='L-BFGS-B',
               options={'maxiter': MAXITER, 'maxcor': 500, 'iprint': 5}, tol=1e-15, 
               callback=callback)

B_dot_n = np.sum(bs.B().reshape((nphi, ntheta, 3)) * s.unitnormal(), axis=2)
print('starting max|B dot n|:', B_dot_n_max)
print('Final max|B dot n|:', np.max(np.abs(B_dot_n)))
#total_curve_length = sum(func.J() for func in Jls)
#print("starting Sum of lengths of base curves:", total_curve_length_start)
#print("Sum of lengths of base curves:", total_curve_length)

plt.plot(call_result)
plt.xlabel('iterations')
plt.ylabel('F(x)')
plt.yscale('log')
plt.grid('minor')
plt.savefig('cost_function.png')
plt.show()

fig = plot(coils + [s] + [surface_coils], engine="plotly", close=True)
fig = plot(coils + [s], engine="plotly", close=True)
fig = plot(coils + [surface_plot_coils], engine="plotly", close=True)

print('currents:')
print([str(coils[i].current.get_value()) for i,_ in enumerate(coils)])

#import pdb; pdb.set_trace()