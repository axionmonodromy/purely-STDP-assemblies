from matplotlib import pyplot as plt
from stdp_assemblies import *
import pars

params = pars.pars_avg_dw

print("Simulating")
dw_numerical = [delta_W(N, params["w_max"], T=params["tend"]) for N in range(4, 17)]
dw_theoretical = [delta_W_th(N, params["w_max"]) for N in np.linspace(4, 17, 100)]

print("Creating plots")

import matplotlib as mpl

mpl.rcParams["text.usetex"] = True
mpl.rcParams.update({"font.size": 20})

fig, ax = plt.subplots(figsize=(7, 6))
ax.plot(range(4, 17), dw_numerical, ".")
ax.plot(np.linspace(4, 17, 100), dw_theoretical)
ax.set_ylim(-1e-6, 6e-6)
plt.xlabel("$N$")
plt.ylabel("$\overline{\Delta W_{ij}} \ (\mathrm{s}^{-1})$")
fig.savefig("results/time_averaged_dw/time_avg_dw.png")
print("Finished")
