import freud
import gsd.hoomd
import matplotlib.pyplot as plt
import numpy as np

bins = 100
k_max = 30
k_min = 3
sfDirect = freud.diffraction.StaticStructureFactorDirect(
    bins=bins, k_max=k_max, k_min=k_min
)
sfDebye = freud.diffraction.StaticStructureFactorDebye(
    num_k_values=bins, k_max=k_max, k_min=k_min
)

with gsd.hoomd.open("data/LJsampletraj.gsd", "r") as traj:
    for frame in traj:
        sfDebye.compute(frame, reset=False)
        sfDirect.compute(frame, reset=False)

plt.plot(sfDebye.k_values, sfDebye.S_k, label="Debye")
plt.plot(sfDirect.bin_centers, sfDirect.S_k, label="Direct")
plt.title("Static Structure Factor")
plt.xlabel("$k$")
plt.ylabel("$S(k)$")
plt.legend()
plt.show()
