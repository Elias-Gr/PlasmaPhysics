""" #%%
from simsopt.geo import CurveXYZFourier

curve = CurveXYZFourier(300,1)
curve.set_dofs([0,0,0,0,1,0,0,0,1])
curve.plot() """

import matplotlib.pyplot as plt
from simsopt._core import load
# replace "NAME_OF_FILE_YOU_DOWNLOADED" with the name you gave the file
res = load(f'serial0021326.json')

coils = res[1]
surface = res[0]


from simsopt.geo import plot
plot(coils + surface, engine="plotly", close=True)

