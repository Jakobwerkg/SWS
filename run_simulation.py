import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D
from datetime import datetime
from config import *
from simulation import rk4_step, upwind_step

# Select time-stepping function
if numerical_method == "RK4":
    step_function = rk4_step
elif numerical_method == "upwind_space_eueler_time":
    step_function = upwind_step

else:
    raise ValueError(f"Unknown method: {numerical_method}")

# Initialize fields
h = np.zeros((nx, ny))
u = np.zeros((nx, ny))
v = np.zeros((nx, ny))

# Gaussian bump in the center
h += bump_amplitude * np.exp(-((X * 1000 - Lx/2)**2 + (Y * 1000 - Ly/2)**2) / (2 * (Lx / bump_width_factor)**2))

# Setup 2 separate figures
fig_h = plt.figure(figsize=(10, 8))
ax_h = fig_h.add_subplot(111, projection='3d')

fig_quiver, ax_quiver = plt.subplots(figsize=(8, 7))

def animate_all(i):
    global h, u, v
    h, u, v = step_function(h, u, v)
    t_min = i * dt / 60

    # --- Height Surface Plot ---
    ax_h.clear()
    surf = ax_h.plot_surface(X, Y, h, cmap='plasma', vmin=-20, vmax=20, rstride=1, cstride=1)
    ax_h.set_title(f"h(x, y) at t = {t_min:.1f} min")
    ax_h.set_xlabel("x (km)")
    ax_h.set_ylabel("y (km)")
    ax_h.set_zlabel("h (m)")
    ax_h.set_zlim(zmin, zmax)

    # --- Quiver Plot of Wind Vectors ---
    ax_quiver.clear()
    skip = 6

    # Plot the background wind speed magnitude (as shading)
    mesh = ax_quiver.pcolormesh(X, Y, np.sqrt(u**2 + v**2)*2, shading='auto', cmap='viridis')

    # Overlay the quiver vectors
    ax_quiver.quiver(X[::skip, ::skip], Y[::skip, ::skip],
                     u[::skip, ::skip], v[::skip, ::skip], scale=40, color='white')

    ax_quiver.set_title(f"Wind Vectors (u, v) at t = {t_min:.1f} min")
    ax_quiver.set_xlabel("x (km)")
    ax_quiver.set_ylabel("y (km)")
    ax_quiver.set_xlim(X.min(), X.max())
    ax_quiver.set_ylim(Y.min(), Y.max())

    # Terminal progress
    percent_done = (i + 1) / nt * 100
    bar = '█' * int(percent_done // 2.5) + '-' * (40 - int(percent_done // 2.5))
    sys.stdout.write(f'\rProgress: |{bar}| {percent_done:6.2f}% ({i + 1}/{nt} iterations)')
    sys.stdout.flush()

    return []

def log(msg):
    now = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(f"[{now}] {msg}")

if __name__ == "__main__":
    start_time = datetime.now()
    log("Simulation started.")
    os.makedirs(output_dir, exist_ok=True)
    log(f"Output directory: {output_dir}")
    log(f"Running {nt} iterations with {numerical_method} method...")

    ani_h = animation.FuncAnimation(fig_h, animate_all, frames=nt, interval=interval, blit=False)
    ani_quiver = animation.FuncAnimation(fig_quiver, animate_all, frames=nt, interval=interval, blit=False)

    # Construct reproducible filenames
    meta_suffix = f"_{numerical_method}_nx{nx}_ny{ny}_dt{dt}_nt{nt}.gif"
    height_filename = os.path.join(output_dir, "height_surface" + meta_suffix)
    wind_filename = os.path.join(output_dir, "wind_vectors" + meta_suffix)

    # Save animations
    ani_h.save(height_filename, writer="pillow", fps=fps)
    ani_quiver.save(wind_filename, writer="pillow", fps=fps)

    log(f"Saved height surface animation to: {height_filename}")
    log(f"Saved wind vector animation to: {wind_filename}")

    end_time = datetime.now()
    log(f"Simulation and animation complete at {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
    log(f"Total simulation and animation time: {end_time - start_time}")
