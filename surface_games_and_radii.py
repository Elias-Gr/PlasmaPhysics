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
surface_coils = SurfaceRZFourier()

MAXITER = 300


surface_coils.set_rc(0,0,1)
surface_coils.set_rc(1,0,0.6)
surface_coils.set_zs(1,0,0.6)

nphi = 32
ntheta = 32
s = SurfaceRZFourier.from_vmec_input('input', range="full torus", nphi=nphi, ntheta=ntheta)

ncoils = 3
base_curves = create_equally_spaced_curves(ncoils, s.nfp, stellsym=True, R0=1, R1=0.7, order=5)
base_currents = [Current(1.0) * 1e5 for i in range(ncoils)]
base_currents[0].fix_all()

coils = coils_via_symmetries(base_curves, base_currents, s.nfp, True)

#plot([surface_coils] + coils +[s], engine="plotly", close=True)

bs = BiotSavart(coils)
bs.set_points(s.gamma().reshape((-1, 3)))

##setting up parameters for penalty function:

CURVATURE_THRESHOLD = 10
DIST_THRESHOLD = 10

Jf = SquaredFlux(s, bs)

Jcc1 = [LpCurveCurvature(c, 1) for c in base_curves]
Jcc = sum(QuadraticPenalty(J, CURVATURE_THRESHOLD, 'max') for J in Jcc1)


### create child class from CurveSurfaceDistance and change penalty to large values of distance

class CurveDistInverse(CurveSurfaceDistance):
    def __init__(self, curves, surface, minimum_distance):
        super().__init__(curves, surface, minimum_distance)

    def compute_candidates(self):
        if self.candidates is None:
            candidates = sopp.get_pointclouds_closer_than_threshold_between_two_collections(
                [c.gamma() for c in self.curves], [self.surface.gamma().reshape((-1, 3))], self.minimum_distance)
            self.candidates = candidates

    def shortest_distance_among_candidates(self):
        self.compute_candidates()
        from scipy.spatial.distance import cdist
        xyz_surf = self.surface.gamma().reshape((-1, 3))
        return np.max([np.max(cdist(self.curves[i].gamma(), xyz_surf)) for i, _ in self.candidates])

    def shortest_distance(self):   #actually the largest but idk where this function is called
        self.compute_candidates()
        if len(self.candidates) > 0:
            return self.shortest_distance_among_candidates()
        from scipy.spatial.distance import cdist
        xyz_surf = self.surface.gamma().reshape((-1, 3))
        return np.max([np.max(cdist(self.curves[i].gamma(), xyz_surf)) for i in range(len(self.curves))])

         
    def J(self):
        """
        This returns the value of the quantity.
        """
        self.compute_candidates()
        res = 0
        gammas = self.surface.gamma().reshape((-1, 3))
        ns = self.surface.normal().reshape((-1, 3))
        for i, _ in self.candidates:
            gammac = self.curves[i].gamma()
            lc = self.curves[i].gammadash()
            res += self.J_jax(gammac, lc, gammas, ns)
        return 1/(res+1e-9)

#import pdb; pdb.set_trace()

Jdist = CurveDistInverse(base_curves, surface_coils, DIST_THRESHOLD)

WEIGHT_DIST = 0.001
WEIGHT_CURVE = 1

JF = Jf + WEIGHT_CURVE * Jcc + WEIGHT_DIST * Jdist


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


res = minimize(fun, dofs, jac=True, method='L-BFGS-B',
               options={'maxiter': MAXITER, 'maxcor': 500, 'iprint': 5}, tol=1e-15)

B_dot_n = np.sum(bs.B().reshape((nphi, ntheta, 3)) * s.unitnormal(), axis=2)
print('starting max|B dot n|:', B_dot_n_max)
print('Final max|B dot n|:', np.max(np.abs(B_dot_n)))
#total_curve_length = sum(func.J() for func in Jls)
#print("starting Sum of lengths of base curves:", total_curve_length_start)
#print("Sum of lengths of base curves:", total_curve_length)


fig = plot(coils + [s] + [surface_coils], engine="plotly", close=True)
fig = plot(coils + [s], engine="plotly", close=True)

#import pdb; pdb.set_trace()