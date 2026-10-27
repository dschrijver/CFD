import h5py
from pathlib import Path
from natsort import natsorted
import matplotlib.pyplot as plt
import numpy as np

data_directory = __file__.rsplit("/", 2)[0] + "/"

data_u_path = Path(data_directory)
data_u_files = data_u_path.glob("u_*.h5")
data_u_files = natsorted(data_u_files, key=lambda x: x.name)
n_files = len(data_u_files)

data_v_path = Path(data_directory)
data_v_files = data_v_path.glob("v_*.h5")
data_v_files = natsorted(data_v_files, key=lambda x: x.name)

with h5py.File(data_u_files[-1], "r") as data_file:
    NX, NY, NZ = data_file["data"].shape
    u = data_file["data"][:,:,0]
    time = data_file["time"][()]

with h5py.File(data_v_files[-1], "r") as data_file:
    NX, NY, NZ = data_file["data"].shape
    v = data_file["data"][:,:,0]
    time = data_file["time"][()]

u_extend = np.zeros((NX + 1, NY))
v_extend = np.zeros((NX, NY + 1))

u_extend[:-1, :] = u
v_extend[:, :-1] = v

u_center = np.empty((NX, NY))
v_center = np.empty((NX, NY))
for i in range(NX):
    for j in range(NY):
        u_center[i,j] = 0.5*(u_extend[i, j] + u_extend[i + 1, j])
        v_center[i,j] = 0.5*(v_extend[i, j] + v_extend[i, j + 1])

dx = 1/NX

x = np.arange(0.5, NX, 1) * dx
y = np.arange(0.5, NY, 1) * dx

velocity = np.sqrt(u_center*u_center + v_center*v_center)
im = plt.imshow(velocity.T, origin="lower", cmap="jet", extent=[0, NX*dx, 0, NY*dx], vmin=0, vmax=1)
plt.streamplot(x, y, u_center.T, v_center.T, color="white", density=1.5, linewidth=1)
plt.colorbar(im, label="velocity")
plt.xlabel("$x$")
plt.ylabel("$y$")
plt.savefig("cavity_streamplot.pdf")
plt.close()

y_exact = [0, 0.0547, 0.0625, 0.0703, 0.1016, 0.1719, 0.2813, 0.4531, 0.5, 0.6172, 0.7344, 0.8516, 0.9531, 0.9609, 0.9688, 0.9766, 1]
u_exact = [0, -0.03717, -0.04192, -0.04775, -0.06434, -0.10150, -0.15662, -0.21090, -0.20581, -0.13641, 0.00332, 0.23151, 0.68717, 0.73722, 0.78871, 0.84123, 1]

x_exact = [0, 0.0625, 0.0703, 0.0781, 0.0938, 0.1563, 0.2266, 0.2344, 0.5, 0.8047, 0.8594, 0.9063, 0.9453, 0.9531, 0.9609, 0.9688, 1]
v_exact = [0, 0.09233, 0.10091, 0.10890, 0.12317, 0.16077, 0.17507, 0.17527, 0.05454, -0.24533, -0.22445, -0.16914, -0.10313, -0.08864, -0.07391, -0.05906, 0]

plt.plot(y, u[NX//2, :], zorder=0, label="CFD")
plt.scatter(y_exact, u_exact, c="red", zorder=1, label="Ghia 1982")
plt.xlabel("$y$")
plt.ylabel("$u$")
plt.suptitle("$u$ profile along vertical in the center of the cavity")
plt.legend()
plt.savefig("cavity_u.pdf")
plt.close()

plt.plot(x, v[:, NY//2], zorder=0, label="CFD")
plt.scatter(x_exact, v_exact, c="red", zorder=1, label="Ghia 1982")
plt.xlabel("$x$")
plt.ylabel("$v$")
plt.suptitle("$v$ profile along horizontal in the center of the cavity")
plt.legend()
plt.savefig("cavity_v.pdf")