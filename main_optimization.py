import matplotlib.pyplot as plt
import numpy as np
import os
import sys
from simsopt._core import load
from simsopt.geo import plot, CurveSurfaceDistance, CurveCurveDistance, LpCurveCurvature
from simsopt.geo import SurfaceRZFourier, create_equally_spaced_curves, curves_to_vtk
from simsopt.field import Current, coils_via_symmetries, BiotSavart
from simsopt.objectives import SquaredFlux, QuadraticPenalty
from scipy.optimize import minimize
import simsoptpp as sopp
from scipy.optimize import OptimizeResult
from parameters import *
from simsopt._core import save
import runPoincare as rp

### change into folder (USE YOUR OWN DIRECTORY)

os.chdir(directory)

if not os.path.exists(foldername):
    os.mkdir(foldername)


#### ---------------- starting conditions and surfaces ------------------ ####

### creating all surfaces to pin the coils on Torus with fixed Radius via CurveSurfaceDistance ####

inner_torus = SurfaceRZFourier()
outer_torus = SurfaceRZFourier()
coil_torus = SurfaceRZFourier()

### change the shape of the surfaces:

inner_torus.set_rc(0,0,R)
inner_torus.set_rc(1,0,r_coil_surface-r_distance_from_surface)
inner_torus.set_zs(1,0,r_coil_surface-r_distance_from_surface)

outer_torus.set_rc(0,0,R)
outer_torus.set_rc(1,0,r_coil_surface+r_distance_from_surface)
outer_torus.set_zs(1,0,r_coil_surface+r_distance_from_surface)

coil_torus.set_rc(0,0,R)
coil_torus.set_rc(1,0,r_coil_surface)
coil_torus.set_zs(1,0,r_coil_surface)

### do another surface not as a torus but a ovaloid on the outside to prevent coils outside the outer torus reach

outer_surface = SurfaceRZFourier()

outer_surface.set_rc(0,0,0)
outer_surface.set_rc(1,0,Radius_outer_surface)
outer_surface.set_zs(1,0,Radius_outer_surface_z)

### importing plasma Surface from vmec

s = SurfaceRZFourier.from_vmec_input('input', range="full torus", nphi=nphi, ntheta=ntheta)

### creating coils as starting condition

if PRESIM_OPT:
    coils_presim = load(f'serial0021326.json')[1]
    #coils_presim = load('test_presim/bs.json').coils
    base_curves = [x.curve for x in coils_presim]
    base_currents = [x.current for x in coils_presim]
else:
    base_curves = create_equally_spaced_curves(ncoils, s.nfp, stellsym=True, R0=1, R1=r_coil_surface, order=fourierordercoils)
    base_currents = [Current(1.0) * CURRENT for i in range(ncoils)]

if FIXEDCURRENT:
    base_currents[0].fix_all()  # has to be used when going for more iterations (~>500)#

coils = coils_via_symmetries(base_curves, base_currents, s.nfp, True)

### plot coils and surface:

#plot([surface_coils] + coils +[s], engine="plotly", close=True, opacity = 0.5)

### creating biotsavard object for flux penalty

bs = BiotSavart(coils)    # calculates magnetic field from coils
bs.set_points(s.gamma().reshape((-1, 3)))



####  ------- cost function ------- ####

Jf = SquaredFlux(s, bs, definition = 'quadratic flux')
Jcc1 = [LpCurveCurvature(c, 1) for c in base_curves]
Jcc = sum(QuadraticPenalty(J, CURVATURE_THRESHOLD, 'max') for J in Jcc1)
Jdist_1 = CurveSurfaceDistance(base_curves, inner_torus, DIST_THRESHOLD)
Jdist_2 = CurveSurfaceDistance(base_curves, outer_torus, DIST_THRESHOLD)
Jdist_3 = CurveSurfaceDistance(base_curves, outer_surface, DIST_THRESHOLD_OUT)
Jccdist = CurveCurveDistance(base_curves, CCDIST_THRESH)



JF = (BDOTN_WEIGHT*Jf + WEIGHT_CURVE * Jcc + WEIGHT_DIST * Jdist_1 + WEIGHT_DIST * Jdist_2 
    + CCDIST_WEIGHT * Jccdist + WEIGHT_DIST_OUT * Jdist_3)



#### create output file

stdoutOrigin=sys.stdout 
sys.stdout = open(foldername+"/log.txt", "w")



#### ----------- minimization ------------ ####


B_dot_n = np.sum(bs.B().reshape((nphi, ntheta, 3)) * s.unitnormal(), axis=2)


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


####  creating function to extract cost function at each iteration step and save into call_result array

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



#### plotting cost function over iteration

plt.plot(call_result)
plt.xlabel('iterations')
plt.ylabel('F(x)')
plt.yscale('log')
plt.grid('minor')
plt.savefig(foldername+'/cost_function.png')



#### plotting coils with different surfaces to inspect complexity and confirm radial pinning


fig = plot(coils + [s], engine="plotly", close=True)
fig.write_html(foldername+'/coils_and_surf.html')
fig = plot(coils + [coil_torus], engine="plotly", close=True)
fig.write_html(foldername+'/coils_fixedradius.html')

#### save coil curves

save(bs,foldername+'/bs.json')

curves = [c.curve for c in coils]
curves_to_vtk(curves, foldername+"/curves_opt")

#### make sure currents don't shrink too much when disabling fix currents

print('currents:')
print([str(coils[i].current.get_value()) for i,_ in enumerate(coils)])

os.chdir(directory)

#### close output file

sys.stdout.close()
sys.stdout=stdoutOrigin

