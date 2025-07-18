"""
plot.py

Plots simulation data in both a static plot and animation to observe the surface
field evolve.

2025.07.18
Anthony Atkinson
"""

import numpy as np
from numpy import ma

from parse_config import parse_config

# very small number
eps = 1e-10

# config file location
cfname = "config.c"
# read config file for params
params = parse_config(cfname, debug=1)

# get params from config file
fname = params.fname
ntheta = params.ntheta
nphi = params.nphi
dth = params.dth
dph = params.dph
nt = params.nt
dt = params.dt
freq = params.freq
nframes = (nt - 1) // freq + 1
tcycle = params.tcycle

data = np.loadtxt(fname)
a = data.reshape(nframes, ntheta, nphi)

theta = np.linspace(-dth, np.pi+dth+eps, num=ntheta, endpoint=True)
theta = theta.reshape(1, ntheta, 1)
sintheta = np.sin(theta).reshape(1, ntheta, 1)
phi = np.linspace(0, 2 * np.pi + 1e-14, num=nphi, endpoint=True)

# flux density
fluxd = a * sintheta * dth * dph
# flux density averaged over space for each timestep
fluxd_mean = np.mean(fluxd, axis=(1,2))

cycleframes = int(tcycle // (dt * freq))
# average over complete stellar cycles
ind_cycle = int((nframes // cycleframes) * cycleframes)
# stellar min is every half cycle
inds_min = np.arange(0, nframes, cycleframes / 2, dtype=np.int64)
# stellar max is every quarter and three-quarter cycle
inds_max = np.arange(cycleframes / 4, nframes, cycleframes / 2, dtype=np.int64)

print(nframes)
print(cycleframes)
print(ind_cycle)
print(inds_min)
print(inds_max)

ra_sol_max = 20             # solar radii
fluxd_sol_max = 3.81        # gauss

fluxd_mean_cycle = np.fabs(np.mean(fluxd_mean[:ind_cycle]))
fluxd_mean_min = np.fabs(np.mean(fluxd_mean[inds_min]))
fluxd_mean_max = np.fabs(np.mean(fluxd_mean[inds_max]))

ra_cycle = ra_sol_max * (fluxd_mean_cycle / fluxd_sol_max)**(-0.16)
ra_min = ra_sol_max * (fluxd_mean_min / fluxd_sol_max)**(-0.16)
ra_max = ra_sol_max * (fluxd_mean_max / fluxd_sol_max)**(-0.16)

print(f"Mean Alfvén Radius (Full Cycle): {ra_cycle:.3f} Rsun")
print(f"Mean Alfvén Radius (Min): {ra_min:.3f} Rsun")
print(f"Mean Alfvén Radius (Max): {ra_max:.3f} Rsun")
