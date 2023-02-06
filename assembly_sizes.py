from stdp_assemblies import *
import pars
import multiprocessing as mp
from scipy import optimize as opt

nsims = 10
wms = np.linspace(0.04, 0.06, 5)
wms_plot = np.linspace(0.04, 0.06, 51)
pars_list = [pars.pars_sizes.copy() for _ in wms]
for i, params in enumerate(pars_list):
    params["w_max"] = wms[i]


def minus_deltaw(N, w_max):
    return -delta_W_th(N, w_max)


max_wm = [
    opt.minimize_scalar(
        minus_deltaw, args=(wm,), bounds=(1, 1 + 1 / wm), method="bounded"
    )["x"]
    for wm in wms_plot
]
zero_wm = [
    opt.root_scalar(minus_deltaw, bracket=(2, 1 + 1 / wm - 1e-5), args=(wm,)).root
    for wm in wms_plot
]


wm_sizes = []
wm_sizes_corr = []
if __name__ == "__main__":
    for params in pars_list:
        pool = mp.Pool(nsims)
        p = [params.copy() for _ in range(nsims)]
        print("Simulating for w_max = {}".format(params["w_max"]))
        n_sim_sizes = pool.map(sizes_sim, p)
        pool.close()
        wm_sizes.append(sum([s[0] for s in n_sim_sizes], []))
        wm_sizes_corr.append(sum([s[1] for s in n_sim_sizes], []))


medians_louvain = [np.median(s) for s in wm_sizes]
errbars_louvain_lower = [np.percentile(s, 25) - np.median(s) for s in wm_sizes]
errbars_louvain_upper = [np.median(s) - np.percentile(s, 75) for s in wm_sizes]
errbars_louvain = np.vstack((errbars_louvain_lower, errbars_louvain_upper))

medians_corr = [np.median(s) for s in wm_sizes_corr]
errbars_corr_lower = [np.percentile(s, 25) - np.median(s) for s in wm_sizes_corr]
errbars_corr_upper = [np.median(s) - np.percentile(s, 75) for s in wm_sizes_corr]
errbars_corr = np.vstack((errbars_corr_lower, errbars_corr_upper))


print("Creating plots")

import matplotlib as mpl

mpl.rcParams["text.usetex"] = True
mpl.rcParams.update({"font.size": 20})


fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111)
ax.errorbar(
    wms,
    medians_louvain,
    yerr=errbars_louvain,
    linestyle="",
    marker=".",
    label="observation",
)
ax.errorbar(
    wms,
    medians_corr,
    yerr=errbars_corr,
    linestyle="",
    marker=".",
    label="corrected observation",
)
ax.plot(wms_plot, zero_wm, label="$N: \ \overline{\Delta W}(N)=0$")
ax.plot(
    wms_plot,
    max_wm,
    label="$\mathrm{argmax}_N \overline{\Delta W}(N)$",
    color="tab:brown",
)
ax.plot(wms_plot, 1 + 1 / wms_plot, label="$1+1/\hat{w}$", color="black")
ax.legend()
plt.xlabel("$\hat w$")
plt.ylabel("$N$")
fig.savefig("results/assembly_sizes/sizes.png")
