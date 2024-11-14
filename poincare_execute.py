import runPoincare as rp
from parameters import *
from simsopt._core import load
from simsopt.geo import SurfaceRZFourier

bs = load(foldername+'/bs.json')

coils = bs.coils
s = SurfaceRZFourier.from_vmec_input('input', range="full torus", nphi=nphi, ntheta=ntheta)

FILEDLINE_PARAMS["EXTCUR"] = [c.current.get_value() for c in coils]

rp.run_poincare_simsopt(coils,'variant',foldername+'/poincare',s,True,FILEDLINE_PARAMS)
rp.plot_poincare(foldername+'/poincare', PLOT_PARAMS, s.nfp, 'variant', save=None)

