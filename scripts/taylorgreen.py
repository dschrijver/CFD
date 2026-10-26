import h5py
from pathlib import Path
from natsort import natsorted
import matplotlib.pyplot as plt
import numpy as np

def analytical_u(x, y, Re, t):
    return np.sin(x)*np.cos(y)*np.exp(-2/Re*t)

data_directory = __file__.rsplit("/", 2)[0] + "/"

data_path = Path(data_directory)
data_files = data_path.glob("u_*.h5")
data_files = natsorted(data_files, key=lambda x: x.name)
n_files = len(data_files)

with h5py.File(data_files[-1], "r") as data_file:
    NX, NY, NZ = data_file["data"].shape
    data = data_file["data"][0,:,0]
    time = data_file["time"][()]

dx = 0.01
Re = 100

y = np.arange(0.5, NY, 1) * dx
LY = NY * dx

plt.plot(data, y)
plt.plot(analytical(y, LY, Re, dpdx), y, ls="--")
plt.xlabel("$u$")
plt.ylabel("$y$")
plt.savefig("poiseuille.pdf")