import numpy as np
from numpy import ma
from matplotlib import pyplot as plt
import matplotlib.animation as anim

ntheta = 128
nphi = 256

nt = 10001
dt = 3e3 / 86400    # days
freq = 100     # steps

field_rad = 6.957e10
field_b0 = 8.5
bmr_b0 = 10.0
bmr_th = (90 - 17.5) / 180 * np.pi
bmr_loc = bmr_th / np.pi * ntheta

fname = "bfld.dat"
data = np.loadtxt(fname)
a = data.reshape(-1, ntheta, nphi)

dth = np.pi / (ntheta - 2)
theta = np.linspace(-dth, np.pi+dth+1e-14, num=ntheta, endpoint=True)
theta = theta.reshape(1, ntheta, 1)
dphi = 2 * np.pi / (nphi - 1)
phi = np.linspace(0, 2 * np.pi + 1e-14, num=nphi, endpoint=True)
dS = dth * dphi * np.abs(np.sin(theta))
field_flux0 = field_b0 * dphi * dth
field_flux0 = field_b0 * dphi * dth * np.sin(dth)

#a = dS * a

nframes = a.shape[0]
ms = 100

a_longavg = ma.array(np.mean(a, axis=2))

fig, ((ax1, ax1l), (ax2, ax2l)) = \
        plt.subplots(ncols=2, nrows=2, figsize=(8, 8), layout="constrained")
im1 = ax1.imshow(a[0], vmin=-field_b0, vmax=field_b0, norm="symlog")
ax1.axhline(y=bmr_loc, color="orange", ls="--", lw=1.0)
ax1.axhline(y=ntheta-bmr_loc, color="orange", ls="--", lw=1.0)
cb1 = plt.colorbar(im1, shrink=1.0, orientation="horizontal")
cb1.ax.set_title("Surface Magnetic Field (G)")
cb1.ax.set_xticks([-field_b0, 0.0, field_b0])
cb1.ax.xaxis.set_major_formatter("{x:.02e}")
ax1.set_xlabel(r"Azimuth $\phi$")
ax1.set_xticks([0, nphi//4, nphi/2, 3*nphi//4, nphi-1])
ax1.set_xticklabels(["0", r"$\pi/2$", r"$\pi$", r"$3\pi/2$", r"$2\pi$"])
ax1.set_yticks([0, ntheta/2, ntheta-1])
ax1.set_yticklabels(["0", r"$\pi/2$", r"$\pi$"])
ax1.set_ylabel(r"Colatitude $\theta$")
ax1.set_title("Initial Condition")

im2 = ax2.imshow(a[-1], vmin=-field_b0, vmax=field_b0, norm="symlog")
ax2.axhline(y=bmr_loc, color="orange", ls="--", lw=1.0)
ax2.axhline(y=ntheta-bmr_loc, color="orange", ls="--", lw=1.0)
cb2 = plt.colorbar(im2, shrink=1.0, orientation="horizontal")
cb2.ax.set_title("Surface Magnetic Field (G)")
cb2.ax.set_xticks([-field_b0, 0.0, field_b0])
cb2.ax.xaxis.set_major_formatter("{x:.02e}")
ax2.set_xlabel(r"Azimuth $\phi$")
ax2.set_xticks([0, nphi//4, nphi/2, 3*nphi//4, nphi-1])
ax2.set_xticklabels(["0", r"$\pi/2$", r"$\pi$", r"$3\pi/2$", r"$2\pi$"])
ax2.set_yticks([0, ntheta/2, ntheta-1])
ax2.set_yticklabels(["0", r"$\pi/2$", r"$\pi$"])
ax2.set_ylabel(r"Colatitude $\theta$")
ax2.set_title(f"Final Magnetic Map ({dt * (nt - 1) :.1f} d)")

imlong = ax1l.imshow(a_longavg.T, vmin=-field_b0, vmax=field_b0, norm="symlog")
cb3 = plt.colorbar(imlong, shrink=1.0, orientation="horizontal")
cb3.ax.set_title("Surface Magnetic Field (G)")
cb3.ax.set_xticks([-field_b0, 0.0, field_b0])
cb3.ax.xaxis.set_major_formatter("{x:.02e}")
ax1l.set_xlabel("Steps")
ax1l.set_xlim(xmin=0, xmax=nframes)
ax1l.set_ylabel(r"Colatitude $\theta$")
ax1l.set_yticks([0, ntheta/2, ntheta-1])
ax1l.set_yticklabels(["0", r"$\pi/2$", r"$\pi$"])
ax1l.set_title("Longitudinally-Averaged\nSurface Field over Time")

ax2l.plot(np.ravel(theta), a_longavg[-1], color="C1", alpha=0.75, label="init")
ax2l.plot(np.ravel(theta), a_longavg[0], ls="--", color="C0", label="final")
ax2l.set_xlim(xmin=0, xmax=np.pi)
ax2l.set_xlabel(r"Colatitude $\theta$")
ax2l.set_xticks([0, np.pi/2, np.pi])
ax2l.set_xticklabels(["0", r"$\pi/2$", r"$\pi$"])
ax2l.set_ylabel("Surface Magnetic Flux Density (G)")
ax2l.set_title("Longitudinally-Averaged\nMagnetic Flux Density")
ax2l.legend()

plt.show()

figa, ((axmap, axlmap), (axtotal, axlong)) = \
        plt.subplots(ncols=2, nrows=2, figsize=(8,8), layout="constrained")

imb = axmap.imshow(a[0], vmin=-field_b0, vmax=field_b0, norm="symlog")
axmap.axhline(y=bmr_loc, color="orange", ls="--", lw=1.0)
axmap.axhline(y=ntheta-bmr_loc, color="orange", ls="--", lw=1.0)
cb1a = plt.colorbar(imb, shrink=1.0, orientation="horizontal")
cb1a.ax.set_title("Surface Magnetic Field (G)")
cb1a.ax.set_xticks([-field_b0, 0.0, field_b0])
cb1a.ax.xaxis.set_major_formatter("{x:.02e}")
axmap.set_xlabel(r"Azimuth $\phi$")
axmap.set_xticks([0, nphi//4, nphi/2, 3*nphi//4, nphi-1])
axmap.set_xticklabels(["0", r"$\pi/2$", r"$\pi$", r"$3\pi/2$", r"$2\pi$"])
axmap.set_yticks([0, ntheta/2, ntheta-1])
axmap.set_yticklabels(["0", r"$\pi/2$", r"$\pi$"])
axmap.set_ylabel(r"Colatitude $\theta$")
axmap.set_title(f"Surface Field")

a_longavg[1:] = ma.masked
imlong = axlmap.imshow(a_longavg.T, vmin=-field_b0, vmax=field_b0, norm="symlog")
cb2a = plt.colorbar(imlong, shrink=1.0, orientation="horizontal")
cb2a.ax.set_title("Surface Magnetic Field (G)")
cb2a.ax.set_xticks([-field_b0, 0.0, field_b0])
cb2a.ax.xaxis.set_major_formatter("{x:.02e}")
axlmap.set_xlabel("Steps")
axlmap.set_xlim(xmin=0, xmax=nframes)
axlmap.set_ylabel(r"Colatitude $\theta$")
axlmap.set_yticks([0, ntheta/2, ntheta-1])
axlmap.set_yticklabels(["0", r"$\pi/2$", r"$\pi$"])
axlmap.set_title("Longitudinally-Averaged\nSurface Field over Time")

times = np.linspace(0, (nt - 1) * dt, nframes, endpoint=True)
fluxd = np.sum(a, axis=(1,2))
linetotal = axtotal.plot(0, fluxd[0])[0]

axtotal.set_xlim(xmin=0, xmax=(nt - 1) * dt)
axtotal.set_xlabel("Time (d)")
axtotal.set_ylim(ymin=np.nanmin(fluxd), ymax=np.nanmax(fluxd))
axtotal.set_ylabel("Signed Flux Density (G)")
axtotal.set_title("Total Signed Magnetic Field (G)")

linelong = axlong.plot(np.ravel(theta), a_longavg[0], color="C1", alpha=0.75)[0]
axlong.plot(np.ravel(theta), a_longavg[0], ls="--", color="C0")[0]
axlong.set_xlim(xmin=0, xmax=np.pi)
axlong.set_xlabel(r"Colatitude $\theta$")
axlong.set_xticks([0, np.pi/2, np.pi])
axlong.set_xticklabels(["0", r"$\pi/2$", r"$\pi$"])
axlong.set_ylabel("Surface Magnetic Flux Density (G)")
axlong.set_title("Longitudinally-Averaged\nMagnetic Flux Density")

figtitle = figa.suptitle(f"Step 0 ({0:.1f} d)")

def update(t):
    imb.set_array(a[t])

    a_longavg[t+1:] = ma.masked
    a_longavg.mask[:t+1] = ma.nomask
    imlong.set_array(a_longavg.T)
    
    figtitle.set_text(f"Surface field evolution step {t} ({dt * t * freq:.1f} d)")
    
    linetotal.set_xdata(times[:t])
    linetotal.set_ydata(fluxd[:t])

    linelong.set_xdata(np.ravel(theta))
    linelong.set_ydata(a_longavg[t])
    return imb, imlong, linetotal, linelong, figtitle

ani = anim.FuncAnimation(fig=figa, func=update, frames=nframes, interval=ms)
plt.show()

ani.save(filename="ftcs-imp-schrijver.gif", writer="pillow")
