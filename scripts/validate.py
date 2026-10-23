import h5py
from pathlib import Path
from natsort import natsorted
import matplotlib.pyplot as plt
import numpy as np

data_directory = __file__.rsplit("/", 2)[0] + "/"

data_path = Path(data_directory)
data_files = data_path.glob("pressure_*.h5")
data_files = natsorted(data_files, key=lambda x: x.name)
n_files = len(data_files)

with h5py.File(data_files[0], "r") as data_file:
    NX, NY, NZ = data_file["data"].shape
    data = data_file["data"][:,:,:]

h = 0.01

j = NY//4
k = NZ//4

data = data[:,j,k]

x = np.arange(0.5, NX, 1) * h
y = (j + 0.5) * h
z = (k + 0.5) * h
LX = NX * h
LY = NY * h
LZ = NZ * h

data_analytical = - (np.sin(2*np.pi*x/LX) * np.sin(2*np.pi*y/LY) * np.sin(2*np.pi*z/LZ)) / ((2*np.pi/LX)**2 + (2*np.pi/LY)**2 + (2*np.pi/LZ)**2)

plt.plot(x, data)
plt.plot(x, data_analytical, ls="--")
plt.show()