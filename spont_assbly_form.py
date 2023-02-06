from matplotlib import pyplot as plt
from stdp_assemblies import *
import pars

print("Initializing")
params = pars.pars_spont
network = hawkes_network(**params)
tend = params["tend"]
w0 = np.copy(network.W)
_, w0sorted = cluster(np.copy(network.W))

print("Simulating network")
network.simulate(10000, warmup=True)
network.simulate(tend)


w = np.copy(network.W)
_, wsorted = cluster(network.W)

print("Creating plots")
fig, axs = plt.subplots(2, 2, figsize=(7, 7))

plotmat(axs[0, 0], w0)
plotmat(axs[1, 0], w)
plotmat(axs[0, 1], w0sorted)
plotmat(axs[1, 1], wsorted)
axs[0, 0].title.set_text("t=0")
axs[0, 1].title.set_text("t=0, sorted")
axs[1, 0].title.set_text("t={}".format(tend))
axs[1, 1].title.set_text("t={}, sorted".format(tend))
fig.savefig("results/spont_assembly_form/spont_assbly.png")
print("Finished")
