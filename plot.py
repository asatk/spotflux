import numpy as np
from matplotlib import pyplot as plt

nt = 10

fname = "bfld.dat"
data = np.loadtxt(fname)
a = data.reshape(nt, 128, 256)

fig, (ax1, ax2) = plt.subplots(ncols=1, nrows=2, sharex=True, sharey=True,
                               figsize=(6, 8))
ax1.imshow(a[0])
ax1.set_xlabel(r"$\phi$")
ax1.set_ylabel(r"$\theta$")
ax1.set_title("Initial Condition")

ax2.imshow(a[1])
ax2.set_xlabel(r"$\phi$")
ax2.set_ylabel(r"$\theta$")
ax2.set_title("Some time later")

fig.tight_layout()

plt.show()
