# Import necessary modules
from simsopt.geo import CurveXYZFourier
from simsopt.field import Current, BiotSavart
from simsopt.mhd import Vmec
from simsopt.objectives import LeastSquaresProblem
from simsopt.solve import Optimizer

# 1. Load the target field from a VMEC file (replace 'wout_filename.nc' with your file path)
vmec = Vmec('input')  # Ensure the VMEC file exists with target field data

# 2. Define the coil shapes using Fourier curves and associate currents with them
# Create a list of CurveXYZFourier objects, each representing a coil with Fourier modes
coils = []
for _ in range(3):
    curve = CurveXYZFourier(quadpoints=30, order=2)  # Adjust order and points as needed
    current = Current(1.0)  # Initial current, can be varied if needed
    coils.append((curve, current))  # Each coil is defined by a curve and current

# 3. Set up Biot-Savart field calculation with combined coils
# BiotSavart expects a list of tuples of (curve, current)
biot_savart = BiotSavart([(curve, current) for curve, current in coils])

# 4. Create the objective function to minimize the error between target and generated field
# Here, we use a LeastSquaresProblem to measure the difference in magnetic fields
target_field = vmec.compute_magnetic_field_on_surface()  # Target field on VMEC surface
objective = LeastSquaresProblem.from_tuples([
    (biot_savart.compute, target_field, 1.0)  # Weight of 1.0 for the field match objective
])

# 5. Fix and free specific DOFs
# Example: Fix the zeroth Fourier mode in X to keep coils anchored
for curve, current in coils:
    curve.fix('xc(0)')  # Fix the x cosine mode of order 0 (can be adjusted as needed)

# 6. Run the optimization
optimizer = Optimizer(objective, method='L-BFGS-B')  # L-BFGS-B is a gradient-based optimizer
result = optimizer.optimize()

# 7. Output results
print("Optimization complete. Final objective value:", result.fun)
for curve, _ in coils:
    curve.plot()  # Visualize each coil's optimized shape
