from stdp_assemblies import *
import pars
import matplotlib as mpl

mpl.rcParams["text.usetex"] = True
mpl.rcParams.update({"font.size": 20})

params = pars.pars_sparse

print("Simulating")
dw_sparse = [
    delta_W_sparse(N, params["w_max"], params["density"], T=params["tend"])
    for N in range(4, 20)
]

print("Creating plots")
fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(range(4, 20), dw_sparse, ".")
ax.set_ylim(-1e-6, 6e-6)
plt.xlabel("$N$")
plt.ylabel("$\overline{\Delta W_{ij}} \ (\mathrm{s}^{-1})$")
fig.savefig("results/sparsity/sparsity.png")
print("Finished")
