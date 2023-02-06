from matplotlib import pyplot as plt
from stdp_assemblies import *
import pars

print("Initializing")
params = pars.pars_learn
network = hawkes_network(**params)


assembly_size = params["assembly_size"]
network.W[:assembly_size, :assembly_size] = params["w_max"]
np.fill_diagonal(network.W, 0.0)
network.W[:, -1] = 0.0
network.W[assembly_size : 2 * assembly_size, -1] = params["w_ext"]
network.W[-1, :] = 0.0
network.f0[-1] = params["rho_ext"]
network.set_f0()
w0 = np.copy(network.W)

print("Warmup network")
network.simulate(10000, warmup=True)


print("External input ON")
network.simulate(params["t_ON"], external=1)
w_input = np.copy(network.W)

print("External input OFF")
network.f0[-1] = 0.0
network.set_f0()
network.simulate(params["t_OFF"], external=1)
w_final = np.copy(network.W)


print("Creating plots")
fig, axs = plt.subplots(1, 3, figsize=(13, 3))

plotmat(axs[0], w0[:-1, :-1])
plotmat(axs[1], w_input[:-1, :-1])
plotmat(axs[2], w_final[:-1, :-1])
axs[0].title.set_text("Initial state")
axs[1].title.set_text("After ext. input")
axs[2].title.set_text("Final state")
fig.savefig("results/learn_assembly/learn_assembly.png")

print("Finished")

