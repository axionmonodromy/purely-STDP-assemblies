from matplotlib import pyplot as plt
from stdp_assemblies import *
import pars

import matplotlib as mpl

mpl.rcParams["text.usetex"] = True
mpl.rcParams.update({"font.size": 20})


print("Initializing")
params = pars.pars_drift2
network = hawkes_network(**params)
network.simulate(10000, warmup=True)


w0 = np.copy(network.W)

print("Assembly drift")
for i in range(4):
    network.f0 = np.random.choice(
        [params["rho"], params["rho_small"]],
        size=params["N"],
        p=[1.0 - params["p_switch"], params["p_switch"]],
    )
    network.set_f0()
    network.simulate(params["t_switch"])


print("Final state")
w1 = np.copy(network.W)


network.f0 = params["rho"] * np.ones(params["N"])
network.set_f0()
network.simulate(params["t_switch"])

w2 = np.copy(network.W)

print("Creating plots")
fig, axs = plt.subplots(1, 3, figsize=(13, 3))

plotmat(axs[0], w0)
plotmat(axs[1], w1)
plotmat(axs[2], w2)
axs[0].title.set_text("Initial state")
axs[1].title.set_text("After drift")
axs[2].title.set_text("After resetting base rates")
fig.savefig("results/assembly_drift2/assembly_drift2.png")

print("Finished")
