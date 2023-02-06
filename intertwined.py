from stdp_assemblies import *
import pars

import matplotlib as mpl
from scipy import stats

mpl.rcParams["text.usetex"] = True
mpl.rcParams.update({"font.size": 20})


print("Initializing")
params = pars.pars_intertwined
network = hawkes_network(**params)
t_sim = params["t_sim"]
n_sims = params["n_sims"]
t_eq = params["t_eq"]
w0 = np.copy(network.W)
w_max = params["w_max"]


mat_list = []

W0 = np.copy(network.W)
W1 = np.copy(network.W) / params["w_max"]
W2 = np.copy(W0)
np.fill_diagonal(W2, 1.0)
sum_conn = np.sum(W0)


network.simulate(10000, warmup=True)
print("Initial simulation")
network.simulate(t_eq)

print("Tracking weight matrix")
for _ in tqdm(range(n_sims)):
    network.simulate(t_sim, tqdmsilence=True)
    mat_list.append(np.copy(network.W))


lost_conn_sums = [np.sum(W0 - W1 * mat) for mat in mat_list]
gained_conn_sums = [
    np.sum((np.ones((params["N"], params["N"])) - W1) * mat) for mat in mat_list
]

lost_conn_mats = [(W0[W0 > 0] - mat[W0 > 0]) for mat in mat_list]
gained_conn_mats = [mat[np.isclose(W2, 0.0)] for mat in mat_list]


gain_mat_corr = [
    [
        stats.pearsonr(gained_conn_mats[j * 100], gained_conn_mats[i])[0]
        for i in range(len(mat_list))
    ]
    for j in range(1, 7)
]
loss_mat_corr = [
    [
        stats.pearsonr(lost_conn_mats[j * 100], lost_conn_mats[i])[0]
        for i in range(len(mat_list))
    ]
    for j in range(1, 7)
]


w = np.copy(network.W)


print("Creating plots")

t_plot = np.arange(len(mat_list)) * t_sim * params["seconds_per_unit_time"]
fn_fig = "errgrowth" + ".png"
fig = plt.figure(figsize=(12, 5))
fig.add_subplot(1, 2, 1)
plt.plot(t_plot, gained_conn_sums / sum_conn)
plt.title("gain")
plt.xlabel("$t(\mathrm{s})$")
fig.add_subplot(1, 2, 2)
plt.plot(t_plot, lost_conn_sums / sum_conn)
plt.title("loss")
plt.xlabel("$t(\mathrm{s})$")
plt.tight_layout(pad=0.8)
plt.savefig("results/intertwined/gained_lost.png")


n_sims = len(mat_list)
fn_fig2 = "errcorr" + ".png"
fig2 = plt.figure(dpi=500, figsize=(8, 12))

ax = fig2.add_subplot(2, 1, 1)
for i in range(6):
    plt.plot(t_plot, loss_mat_corr[i])
plt.title("Correlations of lost connections")
plt.ylabel("Pearson corr.")
plt.xlabel("$t(\mathrm{s})$")
plt.plot(t_plot, np.zeros(n_sims), c="Black")
plt.ylim((-0.1, 0.8))

ax = fig2.add_subplot(2, 1, 2)
for i in range(6):
    plt.plot(t_plot, gain_mat_corr[i])
plt.title("Correlation of gained connections")
plt.ylabel("Pearson corr.")
plt.xlabel("$t(\mathrm{s})$")
plt.plot(t_plot, np.zeros(n_sims), c="Black")

plt.ylim((-0.1, 0.8))
plt.tight_layout(pad=0.5)

plt.savefig("results/intertwined/correlations.png")

fig3 = plt.figure(dpi=500, figsize=(5, 5))
plt.imshow(network.W, cmap="Blues")
plt.savefig("results/intertwined/final_state.png")


fig4 = plt.figure(dpi=500, figsize=(5, 5))
plt.imshow(W0, cmap="Blues")
plt.savefig("results/intertwined/initial_state.png")
