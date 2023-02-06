from stdp_assemblies import *
import pars

import matplotlib as mpl

mpl.rcParams["text.usetex"] = True
mpl.rcParams.update({"font.size": 20})


print("Initializing")
params = pars.pars_overlap
network = hawkes_network(**params)
network.simulate(10000, warmup=True)

w0 = np.copy(network.W)

network.f0[10] = params["rho_small"]
network.set_f0()

print("Simulating network")
network.simulate(params["t_end"])

w1 = np.copy(network.W)

print("Creating plots")
fig, axs = plt.subplots(1, 2, figsize=(9, 3))

plotmat(axs[0], w0)
plotmat(axs[1], w1)
axs[0].title.set_text("Initial state")
axs[1].title.set_text("Final state")
fig.savefig("results/overlap/overlap.png")

print("Finished")
