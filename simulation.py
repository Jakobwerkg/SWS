import numpy as np
from config import nx, ny, dx, dy, g, H, dt

def periodic(arr):
    """Apply periodic BCs by wrapping."""
    return np.pad(arr, ((1,1), (1,1)), mode='wrap')

def ddx(f):
    f_p = periodic(f)
    return (f_p[2:, 1:-1] - f_p[0:-2, 1:-1]) / (2 * dx)

def ddy(f):
    f_p = periodic(f)
    return (f_p[1:-1, 2:] - f_p[1:-1, 0:-2]) / (2 * dy)

def dh_dt(h, u, v):
    return -H * (ddx(u) + ddy(v))

def du_dt(h, u, v):
    return -g * ddx(h)

def dv_dt(h, u, v):
    return -g * ddy(h)

def rk4_step(h, u, v):
    k1_h = dh_dt(h, u, v)
    k1_u = du_dt(h, u, v)
    k1_v = dv_dt(h, u, v)

    k2_h = dh_dt(h + 0.5*dt*k1_h, u + 0.5*dt*k1_u, v + 0.5*dt*k1_v)
    k2_u = du_dt(h + 0.5*dt*k1_h, u + 0.5*dt*k1_u, v + 0.5*dt*k1_v)
    k2_v = dv_dt(h + 0.5*dt*k1_h, u + 0.5*dt*k1_u, v + 0.5*dt*k1_v)

    k3_h = dh_dt(h + 0.5*dt*k2_h, u + 0.5*dt*k2_u, v + 0.5*dt*k2_v)
    k3_u = du_dt(h + 0.5*dt*k2_h, u + 0.5*dt*k2_u, v + 0.5*dt*k2_v)
    k3_v = dv_dt(h + 0.5*dt*k2_h, u + 0.5*dt*k2_u, v + 0.5*dt*k2_v)

    k4_h = dh_dt(h + dt*k3_h, u + dt*k3_u, v + dt*k3_v)
    k4_u = du_dt(h + dt*k3_h, u + dt*k3_u, v + dt*k3_v)
    k4_v = dv_dt(h + dt*k3_h, u + dt*k3_u, v + dt*k3_v)

    h_new = h + (dt/6) * (k1_h + 2*k2_h + 2*k3_h + k4_h)
    u_new = u + (dt/6) * (k1_u + 2*k2_u + 2*k3_u + k4_u)
    v_new = v + (dt/6) * (k1_v + 2*k2_v + 2*k3_v + k4_v)

    return h_new, u_new, v_new

# ========== UPWIND SCHEME (Euler in time, upwind in space) ==========

def upwind_x(f, u):
    f_p = periodic(f)
    df = np.zeros_like(f)
    mask = u >= 0
    df[mask] = (f_p[1:-1, 1:-1][mask] - f_p[0:-2, 1:-1][mask]) / dx
    df[~mask] = (f_p[2:, 1:-1][~mask] - f_p[1:-1, 1:-1][~mask]) / dx
    return df

def upwind_y(f, v):
    f_p = periodic(f)
    df = np.zeros_like(f)
    mask = v >= 0
    df[mask] = (f_p[1:-1, 1:-1][mask] - f_p[1:-1, 0:-2][mask]) / dy
    df[~mask] = (f_p[1:-1, 2:][~mask] - f_p[1:-1, 1:-1][~mask]) / dy
    return df

def dh_dt_upwind(h, u, v):
    return -H * (upwind_x(u, u) + upwind_y(v, v))

def du_dt_upwind(h, u, v):
    return -g * upwind_x(h, u)

def dv_dt_upwind(h, u, v):
    return -g * upwind_y(h, v)

def upwind_step(h, u, v):
    h_new = h + dt * dh_dt_upwind(h, u, v)
    u_new = u + dt * du_dt_upwind(h, u, v)
    v_new = v + dt * dv_dt_upwind(h, u, v)
    return h_new, u_new, v_new
