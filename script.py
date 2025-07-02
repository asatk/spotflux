"""
script.py

A skeleton script for manipulating data from the C simulation in Python.

2025.07.01
Anthony Atkinson
"""

import numpy as np
from matplotlib import pyplot as plt

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
freq = params.freq

# number of frames saved
nframes = (nt - 1) // freq + 1

# load surface field data
data = np.loadtxt(fname)
# re-order raw data to be TIME x THETA x PHI
a = data.reshape(nframes, ntheta, nphi)

# theta values at each grid point (from -dt to pi + dt inclusive)
theta_ax = np.linspace(-dth, np.pi+dth+eps, num=ntheta, endpoint=True)
# reshape to multiply with surface field array
theta = theta_ax.reshape(1, ntheta, 1)

# calculate the mean field of the star over all times and angles
mean_field = np.sum(a * np.sin(theta) * dth * dph) / 4 / np.pi

print(f"Mean surface magnetic field: {mean_field:.4e} G")

# min and max field values of initial condition (for consistent colorbars)
a0min = np.min(a[0])
a0max = np.max(a[0])
# display the first frame (initial condition)
fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, sharex=True, sharey=True)
im1 = ax1.imshow(a[0], vmin=a0min, vmax=a0max)
cb1 = plt.colorbar(im1)
ax1.set_title("Initial Condition")
cb1.ax.set_title("Surface Field (G)")

# display the final frame
im2 = ax2.imshow(a[-1], vmin=a0min, vmax=a0max)
cb2 = plt.colorbar(im2)
ax2.set_title("Final State")
cb2.ax.set_title("Surface Field (G)")

fig.tight_layout()

plt.show()
