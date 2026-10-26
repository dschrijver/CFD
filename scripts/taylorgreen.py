import h5py
from pathlib import Path
from natsort import natsorted
import matplotlib.pyplot as plt
import numpy as np

def analytical_u(x, y, LX, LY, Re, t):
    return np.sin(2*np.pi *x / LX)*np.cos(2 * np.pi * y / LY)*np.exp(-2/Re*(2*np.pi/LX)**2*t)

data_directory = __file__.rsplit("/", 2)[0] + "/"

data_path = Path(data_directory)
data_files = data_path.glob("u_*.h5")
data_files = natsorted(data_files, key=lambda x: x.name)
n_files = len(data_files)

with h5py.File(data_files[-1], "r") as data_file:
    NX, NY, NZ = data_file["data"].shape
    data = data_file["data"][:,:,0]
    time = data_file["time"][()]

dx = 0.01
Re = 100

data_analytical = np.empty(data.shape)

x_list = np.arange(0.5, NX, 1) * dx
y_list = np.arange(0.5, NY, 1) * dx
print(x_list)
for i in range(NX):
    for j in range(NY):
        data_analytical[i, j] = analytical_u(x_list[i], y_list[j], 1, 1, Re, time)

print(np.max(data_analytical), np.max(data))