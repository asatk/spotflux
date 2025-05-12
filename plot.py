import numpy as np
from matplotlib import pyplot as plt

ntheta = 128
nphi = 256

nsteps = 1000000
dt = 3.15e7 / 500000 / 86400    # days
freq = 5000     # steps

field_b0 = 8.5e-4
bmr_b0 = 10.0e-4
bmr_th = np.pi / 3
bmr_loc = bmr_th / np.pi * ntheta

fname = "bfld.dat"
data = np.loadtxt(fname)
a = data.reshape(-1, ntheta, nphi)

nt = a.shape[0]

fig, (ax1, ax2) = plt.subplots(ncols=1, nrows=2, sharex=True, sharey=True,
                               figsize=(6, 8))
im1 = ax1.imshow(a[0], vmin=-field_b0, vmax=field_b0, norm="symlog")
cb1 = plt.colorbar(im1, shrink=0.75, orientation="horizontal", pad=0.2)
cb1.ax.set_title("Surface Magnetic Field Strength (T)")
cb1.ax.set_xticks([-field_b0, 0.0, field_b0])
ax1.set_xlabel(r"Azimuth $\phi$")
ax1.set_xticks([0, nphi//4, nphi/2, 3*nphi//4, nphi-1])
ax1.set_xticklabels(["0", r"$\pi/2$", r"$\pi$", r"$3\pi/2$", r"$2\pi$"])
ax1.set_yticks([0, ntheta/2, ntheta-1])
ax1.set_yticklabels(["0", r"$\pi/2$", r"$\pi$"])
ax1.set_ylabel(r"Colatitude $\theta$")
ax1.set_title("Initial Condition")

im2 = ax2.imshow(a[-1], vmin=-field_b0, vmax=field_b0, norm="symlog")
ax2.axhline(y=bmr_loc, color="red", ls="--", lw=1.0)
cb2 = plt.colorbar(im2, shrink=0.75, orientation="horizontal", pad=0.2)
cb2.ax.set_title("Surface Magnetic Field Strength (T)")
cb2.ax.set_xticks([-field_b0, 0.0, field_b0])
ax2.set_xlabel(r"Azimuth $\phi$")
ax2.set_xticks([0, nphi//4, nphi/2, 3*nphi//4, nphi-1])
ax2.set_xticklabels(["0", r"$\pi/2$", r"$\pi$", r"$3\pi/2$", r"$2\pi$"])
ax2.set_yticks([0, ntheta/2, ntheta-1])
ax2.set_yticklabels(["0", r"$\pi/2$", r"$\pi$"])
ax2.set_ylabel(r"Colatitude $\theta$")
ax2.set_title(f"Final Bipolar Magnetic Region ({dt * (nt - 1) * freq:.1f} d)")

fig.tight_layout()

plt.show()

plt.ion()
plt.figure(1, figsize=(6,4))
for t in range(nt):
    plt.clf()
    im = plt.imshow(a[t], vmin=-field_b0, vmax=field_b0, norm="symlog")
    plt.axhline(y=bmr_loc, color="red", ls="--", lw=1.0)
    cb = plt.colorbar(im, shrink=0.75, orientation="horizontal", pad=0.25)
    cb.ax.set_title("Surface Magnetic Field Strength (T)")
    cb.ax.set_xticks([-field_b0, 0.0, field_b0])
    plt.xlabel(r"Azimuth $\phi$")
    plt.ylabel(r"Colatitude $\theta$")
    plt.title(f"Surface field evolution step {t} ({dt * t * freq:.1f} d)")
    plt.xticks([0, nphi//4, nphi/2, 3*nphi//4, nphi-1],
               labels=["0", r"$\pi/2$", r"$\pi$", r"$3\pi/2$", r"$2\pi$"])
    plt.yticks([0, ntheta/2, ntheta-1],
               labels=["0", r"$\pi/2$", r"$\pi$"])
    plt.draw()
    plt.savefig(f"frames/frame{t:02d}.png")

    plt.pause(0.02)
    plt.show()

plt.ioff()
plt.show()
