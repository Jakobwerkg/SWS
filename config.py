import numpy as np
import sys

# ----------- Choose numerical method -----------
numerical_method = "upwind_space_euler_time"  # Options: "RK4", upwind_space_euler_time

valid_methods = ["RK4", "upwind_space_euler_time"]
if numerical_method not in valid_methods:
    raise ValueError(f"Invalid method: {numerical_method}. Choose one of {valid_methods}")

# ----------- Grid and physical setup -----------
nx = ny =  100
Lx = Ly = 1e5   # meters
dx, dy = Lx / nx, Ly / ny
x = np.linspace(0, Lx, nx)
y = np.linspace(0, Ly, ny)
X, Y = np.meshgrid(x / 1000, y / 1000, indexing='ij')  # For plotting in km

# ----------- Physical constants -----------
g = 9.81       # m/s²
H = 1000       # m
c = np.sqrt(g * H)

# ----------- Time stepping -----------
dt = 0.5        # seconds
nt = 5000 # number of time steps

# ----------- CFL condition check -----------
CFL_actual = dt * c / min(dx, dy)
CFL_max = 0.5  # Maximum CFL number for stability

if CFL_actual > CFL_max:
    print(f"\n[WARNING] Your chosen dt = {dt:.3f}s is too large!")
    print(f"           It violates the CFL condition:")
    print(f"           CFL = {CFL_actual:.3f} > CFL_max = {CFL_max}")
    print("           Reduce dt or increase spatial resolution (dx, dy).\n")
    sys.exit(1)
else:
    print(f"[INFO] CFL condition satisfied: CFL = {CFL_actual:.3f} <= {CFL_max}")

# ----------- Initial condition parameters -----------
bump_amplitude = 40  # meters
bump_width_factor = 20  # domain size / bump width

# ----------- Animation settings -----------
zmin, zmax = -50, 50
fps = 28
interval = 50  # ms between frames

# ----------- Output path and file name -----------
output_dir = "out"
output_filename = (
    f"{output_dir}/shallow_water_{numerical_method}_nt{nt}_dt{dt:.0f}_"
    f"{nx}x{ny}.gif"
)
