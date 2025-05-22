import numpy as np
from matplotlib import pyplot as plt

ntheta = 128
nphi = 256

nsteps = 100000
dt = 3e1 / 86400    # days
freq = 10000     # steps

field_rad = 6.957e10
field_b0 = 8.5 * 1000
bmr_b0 = 10.0
bmr_th = np.pi / 3
bmr_loc = bmr_th / np.pi * ntheta

fname = "bfld.dat"
data = np.loadtxt(fname)
a = data.reshape(-1, ntheta, nphi)

dth = np.pi / (ntheta - 2)
dphi = 2 * np.pi / (nphi - 1)
theta = np.linspace(-dth, np.pi+dth+1e-14, num=ntheta)
theta = theta.reshape(1, ntheta, 1)
dS = dth * dphi * np.abs(np.sin(theta))
field_flux0 = field_b0 * dphi * dth
field_flux0 = field_b0 * dphi * dth * np.sin(dth)

#a = dS * a

nt = a.shape[0]

fig, (ax1, ax2) = plt.subplots(ncols=1, nrows=2, sharex=True, sharey=True,
                               figsize=(6, 8))
im1 = ax1.imshow(a[0], vmin=-field_b0, vmax=field_b0, norm="symlog")
#im1 = ax1.imshow(a[0], vmin=-field_flux0, vmax=field_flux0, norm="symlog")
#im1 = ax1.imshow(a[0], vmin=-field_flux0, vmax=field_flux0, norm="linear")
cb1 = plt.colorbar(im1, shrink=0.75, orientation="horizontal", pad=0.2)
cb1.ax.set_title("Surface Magnetic Field (G)")
#cb1.ax.set_title("Surface Magnetic Flux (Mx)")
cb1.ax.set_xticks([-field_b0, 0.0, field_b0])
#cb1.ax.set_xticks([-field_flux0, 0.0, field_flux0])
cb1.ax.xaxis.set_major_formatter("{x:.02e}")
ax1.set_xlabel(r"Azimuth $\phi$")
ax1.set_xticks([0, nphi//4, nphi/2, 3*nphi//4, nphi-1])
ax1.set_xticklabels(["0", r"$\pi/2$", r"$\pi$", r"$3\pi/2$", r"$2\pi$"])
ax1.set_yticks([0, ntheta/2, ntheta-1])
ax1.set_yticklabels(["0", r"$\pi/2$", r"$\pi$"])
ax1.set_ylabel(r"Colatitude $\theta$")
ax1.set_title("Initial Condition")

im2 = ax2.imshow(a[-1], vmin=-field_b0, vmax=field_b0, norm="symlog")
#im2 = ax2.imshow(a[-1], vmin=-field_flux0, vmax=field_flux0, norm="symlog")
ax2.axhline(y=bmr_loc, color="red", ls="--", lw=1.0)
cb2 = plt.colorbar(im2, shrink=0.75, orientation="horizontal", pad=0.2)
cb2.ax.set_title("Surface Magnetic Field (G)")
#cb2.ax.set_title("Surface Magnetic Flux (Mx)")
cb2.ax.set_xticks([-field_b0, 0.0, field_b0])
#cb2.ax.set_xticks([-field_flux0, 0.0, field_flux0])
cb2.ax.xaxis.set_major_formatter("{x:.02e}")
ax2.set_xlabel(r"Azimuth $\phi$")
ax2.set_xticks([0, nphi//4, nphi/2, 3*nphi//4, nphi-1])
ax2.set_xticklabels(["0", r"$\pi/2$", r"$\pi$", r"$3\pi/2$", r"$2\pi$"])
ax2.set_yticks([0, ntheta/2, ntheta-1])
ax2.set_yticklabels(["0", r"$\pi/2$", r"$\pi$"])
ax2.set_ylabel(r"Colatitude $\theta$")
ax2.set_title(f"Final Magnetic Map ({dt * (nt - 1) * freq:.1f} d)")

fig.tight_layout()

plt.show()

r"""
figt, axt = plt.subplots()
imt = axt.imshow(a[49], vmin=-field_flux0, vmax=field_flux0, norm="symlog")
axt.axhline(y=bmr_loc, color="red", ls="--", lw=1.0)
cbt = plt.colorbar(imt, shrink=0.75, orientation="horizontal", pad=0.2)
cbt.ax.set_title("Surface Magnetic Flux (Mx)")
cbt.ax.set_xticks([-field_flux0, 0.0, field_flux0])
axt.set_xlabel(r"Azimuth $\phi$")
axt.set_xticks([0, nphi//4, nphi/2, 3*nphi//4, nphi-1])
axt.set_xticklabels(["0", r"$\pi/2$", r"$\pi$", r"$3\pi/2$", r"$2\pi$"])
axt.set_yticks([0, ntheta/2, ntheta-1])
axt.set_yticklabels(["0", r"$\pi/2$", r"$\pi$"])
axt.set_ylabel(r"Colatitude $\theta$")
axt.set_title(f"Frame 49")
figt.tight_layout()
plt.show()
"""

plt.ion()
plt.figure(1, figsize=(6,4))
for t in range(nt):
    plt.clf()
    im = plt.imshow(a[t], vmin=-field_b0, vmax=field_b0, norm="symlog")
    #im = plt.imshow(a[t], vmin=-field_flux0, vmax=field_flux0, norm="symlog")
    plt.axhline(y=bmr_loc, color="red", ls="--", lw=1.0)
    cb = plt.colorbar(im, shrink=0.75, orientation="horizontal", pad=0.25)
    cb.ax.set_title("Surface Magnetic Field (G)")
    #cb.ax.set_title("Surface Magnetic Flux (Mx)")
    cb.ax.set_xticks([-field_b0, 0.0, field_b0])
    #cb.ax.set_xticks([-field_flux0, 0.0, field_flux0])
    cb.ax.xaxis.set_major_formatter("{x:.02e}")
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
