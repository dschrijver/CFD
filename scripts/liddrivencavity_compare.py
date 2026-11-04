import h5py
from pathlib import Path
from natsort import natsorted
import matplotlib.pyplot as plt
import numpy as np

data_directory = __file__.rsplit("/", 2)[0] + "/"

# N_list = [8, 16, 32, 64]
N_list = [8, 16, 32, 64, 128]
markers = ["o", "s", "^", "D", "a"]

x_exact = [0, 0.0625, 0.0703, 0.0781, 0.0938, 0.1563, 0.2266, 0.2344, 0.5, 0.8047, 0.8594, 0.9063, 0.9453, 0.9531, 0.9609, 0.9688, 1]
y_exact = [0, 0.0547, 0.0625, 0.0703, 0.1016, 0.1719, 0.2813, 0.4531, 0.5, 0.6172, 0.7344, 0.8516, 0.9531, 0.9609, 0.9688, 0.9766, 1]

# # Re = 100
# u_exact = [0, -0.03717, -0.04192, -0.04775, -0.06434, -0.10150, -0.15662, -0.21090, -0.20581, -0.13641, 0.00332, 0.23151, 0.68717, 0.73722, 0.78871, 0.84123, 1]
# v_exact = [0, 0.09233, 0.10091, 0.10890, 0.12317, 0.16077, 0.17507, 0.17527, 0.05454, -0.24533, -0.22445, -0.16914, -0.10313, -0.08864, -0.07391, -0.05906, 0]

# Re = 1000
u_exact = [0, -0.18109, -0.20196, -0.22220, -0.29730, -0.38289, -0.27805, -0.10648, -0.06080, 0.05702, 0.18719, 0.33304, 0.46604, 0.51117, 0.57492, 0.65928, 1]
v_exact = [0, 0.27485, 0.29012, 0.30353, 0.32627, 0.37095, 0.33075, 0.32235, 0.02526, -0.31966, -0.42665, -0.51550, -0.39188, -0.33714, -0.27669, -0.21388, 0]

for n in range(len(N_list)):
    N = N_list[n]
    marker = markers[n]
    with h5py.File(data_directory + "u_100_%d.h5"%N, "r") as data_file:
        NX, NY, NZ = data_file["data"].shape
        u = data_file["data"][:,:,0]

    dx = 1/NX
    y = np.arange(0.5, NY, 1) * dx

    plt.plot(y, u[NX//2, :], zorder=0, label="$N = %d$"%N)
    if N < 32:
        plt.scatter(y, u[NX//2, :], zorder=0, marker=marker)

plt.scatter(y_exact, u_exact, marker='x', zorder=1, color="black", label="Ghia 1982")
plt.xlabel("$y$")
plt.ylabel("$u$")
plt.legend()
plt.savefig("cavity_compare_u.pdf")
plt.close()


for n in range(len(N_list)):
    N = N_list[n]
    marker = markers[n]
    with h5py.File(data_directory + "v_100_" + str(N) + ".h5", "r") as data_file:
        NX, NY, NZ = data_file["data"].shape
        v = data_file["data"][:,:,0]

    dx = 1/NX
    x = np.arange(0.5, NX, 1) * dx
    
    plt.plot(x, v[:, NY//2], zorder=0, label="$N = %d$"%N)
    if N < 32:
        plt.scatter(x, v[:, NY//2], zorder=0, marker=marker)

plt.scatter(x_exact, v_exact, marker='x', zorder=1, color="black", label="Ghia 1982")
plt.xlabel("$x$")
plt.ylabel("$v$")
plt.suptitle("$v$ profile along horizontal in the center of the cavity")
plt.legend()
plt.savefig("cavity_compare_v.pdf")
plt.close()