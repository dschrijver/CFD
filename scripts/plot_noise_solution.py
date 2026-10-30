import h5py
from pathlib import Path
from natsort import natsorted
import matplotlib.pyplot as plt

data_directory = __file__.rsplit("/", 2)[0] + "/test_poisson/"

data_path = Path(data_directory)
data_files = data_path.glob("pressure_*.h5")
data_files = natsorted(data_files, key=lambda x: x.name)
n_files = len(data_files)

with h5py.File(data_files[0], "r") as data_file:
    NX, NY, NZ = data_file["data"].shape
    k = NZ//4
    data = data_file["data"][:,:,k]

h = 0.01
z = (k + 0.5) * h
LX = NX * h
LY = NY * h
LZ = NZ * h

plt.imshow(data.T, origin="lower", extent=[0, LX, 0, LY])
plt.colorbar(label="pressure")
plt.xlabel("$x$")
plt.ylabel("$y$")
plt.suptitle("$z = %.3f$"%z)
plt.savefig("noise_solution.pdf")