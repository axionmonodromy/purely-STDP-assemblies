from matplotlib import pyplot as plt
from stdp_assemblies import *
import pars

print("Initializing")
params = pars.pars_drift
network = hawkes_network(**params)
T1 = params["T1"]
T2 = params["T2"]
w0 = np.copy(network.W)
# _, w0sorted = cluster(np.copy(network.W))

print("First part of simulation")
network.simulate(10000, warmup=True)
network.simulate(T1)
w1 = np.copy(network.W)
L, w1sorted = cluster(network.W)


print("Second part of simulation")
network.simulate(T2)
w2 = sort_matrix(network.W, L)
_, w2sorted = cluster(network.W)


print("Creating plots")

T1 = T1 * pars.seconds_per_unit_time
T2 = T2 * pars.seconds_per_unit_time
fig, axs = plt.subplots(3, 3, figsize=(9, 8))
plotmat(axs[0, 0], w0)
plotmat(axs[1, 0], w1)
plotmat(axs[1, 1], w1sorted)
plotmat(axs[2, 1], w2)
plotmat(axs[2, 2], w2sorted)

axs[0, 0].title.set_text("t=0s")
axs[0, 1].title.set_text("sorted at t={:.0e}s".format(T1))
axs[0, 2].title.set_text("sorted at t={:.0e}s".format(T1 + T2))

axs[0, 0].set_ylabel("t=0s", rotation=0, size="large", labelpad=60)
axs[1, 0].set_ylabel("t={:.0e}s".format(T1), rotation=0, size="large", labelpad=50)
axs[2, 0].set_ylabel("t={:.0e}s".format(T1 + T2), rotation=0, size="large", labelpad=44)

axs[0, 1].axis("off")
axs[0, 2].axis("off")
axs[1, 2].axis("off")
axs[2, 0].spines["top"].set_visible(False)
axs[2, 0].spines["right"].set_visible(False)
axs[2, 0].spines["bottom"].set_visible(False)
axs[2, 0].spines["left"].set_visible(False)
axs[2, 0].set_xticks([])
axs[2, 0].set_yticks([])

for axr in axs:
    for ax in axr:
        ax.set_xticks([])
        ax.set_yticks([])

fig.tight_layout()
fig.savefig("results/assembly_drift/assembly_drift.png")
print("Finished")
