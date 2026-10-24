import h5py
from pathlib import Path
from natsort import natsorted
import matplotlib.pyplot as plt
import numpy as np

def analytical(x, freq):
    return (np.sin(freq*2*np.pi*x/LX) * np.sin(freq*2*np.pi*y/LY) * np.sin(freq*2*np.pi*z/LZ)) / ((freq*2*np.pi/LX)**2 + (freq*2*np.pi/LY)**2 + (freq*2*np.pi/LZ)**2)

data_directory = __file__.rsplit("/", 2)[0] + "/"

data_path = Path(data_directory)
data_files = data_path.glob("pressure_*.h5")
data_files = natsorted(data_files, key=lambda x: x.name)
n_files = len(data_files)

with h5py.File(data_files[0], "r") as data_file:
    NX, NY, NZ = data_file["data"].shape
    data = data_file["data"][:,:,:]

h = 0.01

j = NY//2
k = NZ//2

data = data[:,j,k]

x = np.arange(0.5, NX, 1) * h
y = (j + 0.5) * h
z = (k + 0.5) * h
LX = NX * h
LY = NY * h
LZ = NZ * h

data_analytical = -analytical(x, 1) - 0.5*analytical(x, 10) - 0.1*analytical(x, 20)

plt.plot(x, data)
plt.plot(x, data_analytical, ls="--")
plt.xlabel("$x$")
plt.ylabel("pressure")
plt.suptitle("$y = %.3f, z = %.3f$"%(y, z))
plt.savefig("validation.pdf")