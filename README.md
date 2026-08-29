# SWS — Shallow Water Simulator

A 2D **shallow water equations** solver in NumPy, with animated output: a 3D
surface of the free-surface height and a quiver map of the velocity field.

> Animations are written to `out/`, which is git-ignored (the GIFs run to
> hundreds of MB). Run the simulation to generate them.

---

## The model

The code integrates the **linearized shallow water equations** on a non-rotating
domain with a flat bottom:

$$
\frac{\partial h}{\partial t} = -H\left(\frac{\partial u}{\partial x} + \frac{\partial v}{\partial y}\right)
\qquad
\frac{\partial u}{\partial t} = -g\frac{\partial h}{\partial x}
\qquad
\frac{\partial v}{\partial t} = -g\frac{\partial h}{\partial y}
$$

* `h` — free-surface displacement (m)
* `u, v` — depth-averaged horizontal velocities (m/s)
* `H` — mean fluid depth (m), `g` — gravity (m/s²)

Waves travel at the gravity-wave speed $c = \sqrt{gH}$ (≈ 99 m/s for the default
`H = 1000 m`).

**Domain:** 100 km × 100 km, **periodic** in both directions — waves leaving one
edge re-enter on the opposite side.

**Initial condition:** fluid at rest with a Gaussian bump of surface height in the
middle of the domain — a "dam break in a bathtub" that radiates concentric waves
outward.

---

## Numerical methods

Two schemes live in [`simulation.py`](simulation.py); pick one in
[`config.py`](config.py) via `numerical_method`:

| `numerical_method` | Space | Time | Notes |
|---|---|---|---|
| `"RK4"` | 2nd-order centred differences | classical 4-stage Runge–Kutta | Low dissipation, keeps wave amplitude. Centred differences are dispersive, so short waves ripple. |
| `"upwind_space_eueler_time"` | 1st-order upwind | forward Euler | Robust and monotone, but strongly diffusive — the bump visibly damps over time. |

Spatial derivatives use periodic padding (`np.pad(..., mode='wrap')`), so no
explicit boundary handling is needed anywhere else.

### Stability

`config.py` checks the CFL condition before anything runs:

$$
\mathrm{CFL} = \frac{c\,\Delta t}{\min(\Delta x, \Delta y)} \le 0.5
$$

If `dt` violates it, the run aborts with an explanatory message instead of
producing a blown-up animation.

---

## Getting started

**Requirements:** Python 3.11+, `numpy`, `matplotlib`, `pillow` (GIF writer).

```bash
pip install numpy matplotlib pillow
```

```bash
python run_simulation.py
```

The run prints a progress bar and timestamped log lines, then writes two GIFs to
`out/`, named after the parameters used:

```
out/height_surface_<method>_nx<nx>_ny<ny>_dt<dt>_nt<nt>.gif
out/wind_vectors_<method>_nx<nx>_ny<ny>_dt<dt>_nt<nt>.gif
```

---

## Configuration

Everything tunable lives in [`config.py`](config.py):

| Setting | Default | Meaning |
|---|---|---|
| `numerical_method` | `"upwind_space_eueler_time"` | `"RK4"` or `"upwind_space_eueler_time"` |
| `nx`, `ny` | `100` | Grid points per direction |
| `Lx`, `Ly` | `1e5` | Domain size (m) |
| `H` | `1000` | Mean depth (m) |
| `g` | `9.81` | Gravity (m/s²) |
| `dt` | `0.5` | Time step (s) |
| `nt` | `5000` | Number of time steps |
| `bump_amplitude` | `40` | Initial bump height (m) |
| `bump_width_factor` | `20` | Bump width = domain size / this |
| `zmin`, `zmax` | `-50, 50` | z-axis limits of the surface plot |
| `fps`, `interval` | `28`, `50` | GIF frame rate / frame delay (ms) |
| `output_dir` | `"out"` | Where GIFs are written |

**Rendering is the bottleneck, not the physics.** Each frame redraws a full 3D
surface plus a quiver plot, so `nt = 5000` at `nx = 100` takes a long time and
produces very large GIFs. For a quick look, try `nt = 200`.

---

## Files

| File | Purpose |
|---|---|
| [`config.py`](config.py) | Parameters, CFL check, output naming |
| [`simulation.py`](simulation.py) | Derivative operators and both time-stepping schemes |
| [`run_simulation.py`](run_simulation.py) | Initial condition, animation, GIF export |
| `out/` | Generated animations (git-ignored) |

---

## Notes

* The method string `upwind_space_eueler_time` has a typo ("eueler"), but it is the
  literal value the code compares against — keep it, or rename it in both
  `config.py` and `run_simulation.py`.

---

<sub>README produced by Claude Opus 5.</sub>
