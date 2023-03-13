import numpy as np
import matplotlib.pyplot as plt
import bct
import itertools
from tqdm import tqdm


def stdp_prop_sym(t, tau_p, tau_d):
    return np.array(
        [np.exp(-t / tau_p), np.exp(-t / tau_p), np.exp(-t / tau_d), np.exp(-t / tau_d)]
    )


def stdp_prop_asym(t, tau_p, tau_d):
    return np.array([np.exp(-t / tau_p), np.exp(-t / tau_d)])


def rand_init_weight_matrix(N, w_max, scale=0.20):
    W = np.random.rand(N, N)
    np.fill_diagonal(W, 0)
    W = scale * w_max * W
    return W


def assembly_weight_matrix(N, num_assemblies, assembly_size, w_max):
    W = np.zeros((N, N))
    for i in range(num_assemblies):
        W[
            assembly_size * i : assembly_size * (i + 1),
            assembly_size * i : assembly_size * (i + 1),
        ] = 1.0
    np.fill_diagonal(W, 0.0)
    return w_max * W


def overlap_matrix(nc, swaps):
    olm = np.tri(nc, nc, -1)

    for _ in range(swaps):
        entries = [tuple(x) for x in np.argwhere(olm)]
        a = np.random.choice(len(entries), 2, replace=False)
        olm[entries[a[0]]] += 1
        olm[entries[a[1]]] -= 1
    return olm


def overlap_matrix_remove(nc, removals):
    olm = np.tri(nc, nc, -1)

    for _ in range(removals):
        entries = [tuple(x) for x in np.argwhere(olm)]
        a = np.random.choice(len(entries))
        olm[entries[a]] = 0
    return olm


def intertwined_weight_matrix(N, w_max, nc, swaps):
    #    olm = overlap_matrix(nc,swaps)
    olm = overlap_matrix_remove(nc, swaps)
    neurons = list(range(N))
    assemblies = []
    ass_vec = [np.zeros(N) for _ in range(nc)]
    col_sums = np.sum(olm, axis=0)

    for k in range(nc):

        assembly = []
        assembly = [neurons.pop(0) for _ in range(int(col_sums[k]))]
        for i in range(k):
            for _ in range(int(olm[k][i])):
                assembly.append(assemblies[i].pop(0))

        for i in assembly:
            ass_vec[k][i] = 1

        assemblies.append(assembly)

    W = np.zeros((N, N))
    for vec in ass_vec:
        W += np.outer(vec, vec)
    W = w_max * W
    np.fill_diagonal(W, 0)
    N_reduced = assembly[-1] + 1
    W_reduced = np.copy(W[:N_reduced, :N_reduced])
    return W_reduced


def CDFm1(f1sum, x):
    # inverse cumulative probability density for the time of the next evoked spike in a Hawkes network
    return -np.log(1 + 1 / f1sum * np.log(1 - x))


def cluster(W, louv=1):
    w = np.copy(W)
    try:
        (L, q) = bct.community_louvain(w, louv)
    except:
        return [1 for _ in range(w.shape[0])], w
    else:
        W_ = sort_matrix(w, L)
    return L, W_


def sort_matrix(W, L):
    inds = np.argsort(L)
    W_ = W[:, inds]
    W_ = W_[inds, :]
    return W_


def assembly_sizes(L):
    return [list(L).count(i) for i in range(1, max(L) + 1)]


def plotmat(ax, mat):
    ax.imshow(mat, cmap="Blues")
    return ax


def delta_W(
    N,
    w_max,
    tau_syn=1.0,
    tau_p=2.5,
    tau_d=5.0,
    amp_p=0.08,
    amp_d=-0.0533,
    rho=0.0015,
    T=20_000_000,
    init_W=None,
    seed=None,
    mu=1,
):
    network = hawkes_network(**locals())
    network.W = w_max * np.ones((N, N))
    np.fill_diagonal(network.W, 0.0)
    network.simulate(10000, warmup=True)
    dW = network.track_dw(T)
    return np.sum(dW) / T / (N * (N - 1))


def delta_W_sparse(
    N,
    w_max,
    density,
    tau_syn=1.0,
    tau_p=2.5,
    tau_d=5.0,
    amp_p=0.08,
    amp_d=-0.0533,
    rho=0.0015,
    T=20_000_000,
    init_W=None,
    seed=None,
    mu=1,
):
    network = hawkes_network(**locals())
    network.W = w_max * np.ones((N, N))
    np.fill_diagonal(network.W, 0.0)
    network.Adja = np.random.choice(2, (N, N), p=[1 - density, density])
    np.fill_diagonal(network.Adja, 0.0)
    network.W = network.Adja * network.W

    delta_Ws = []
    network.simulate(10000, warmup=True)
    for _ in range(10):
        W_track = network.track_dw(T)
        delta_Ws.append(np.sum(W_track) / T / (N * (N - 1)))
    return np.mean(delta_Ws), np.std(delta_Ws)


def delta_W_th(
    N, w_max, tau_syn=1.0, tau_p=2.5, tau_d=5.0, amp_p=0.08, amp_d=-0.0533, rho=0.0015
):
    rate_term = (
        2
        * rho
        * rho
        * (amp_p * tau_p + amp_d * tau_d)
        / (1 - (N - 1) * w_max)
        / (1 - (N - 1) * w_max)
    )
    interaction_term = w_max * rho / (1 + w_max) / (1 - (N - 1) * w_max) ** 2 * (
        amp_p
        * tau_p
        * (
            tau_p * (2 + (2 - N) * w_max)
            + tau_syn * (2 - (N - 2) * w_max - (N - 1) * w_max * w_max)
        )
    ) / (tau_syn + (1 + w_max) * tau_p) / (
        tau_syn + tau_p * (1 - (N - 1) * w_max)
    ) + w_max * rho / (
        1 + w_max
    ) / (
        1 - (N - 1) * w_max
    ) ** 2 * (
        amp_d
        * tau_d
        * (
            tau_d * (2 + (2 - N) * w_max)
            + tau_syn * (2 - (N - 2) * w_max - (N - 1) * w_max * w_max)
        )
    ) / (
        tau_syn + (1 + w_max) * tau_d
    ) / (
        tau_syn + tau_d * (1 - (N - 1) * w_max)
    )
    return rate_term + interaction_term


def delta_W_th_nozero(
    N, w_max, tau_syn=1.0, tau_p=2.5, tau_d=5.0, amp_p=0.08, amp_d=-0.0533, rho=0.0015
):
    interaction_term = w_max * rho / (1 + w_max) / (1 - (N - 1) * w_max) ** 2 * (
        amp_p
        * tau_p
        * (
            tau_p * (2 + (2 - N) * w_max)
            + tau_syn * (2 - (N - 2) * w_max - (N - 1) * w_max * w_max)
        )
    ) / (tau_syn + (1 + w_max) * tau_p) / (
        tau_syn + tau_p * (1 - (N - 1) * w_max)
    ) + w_max * rho / (
        1 + w_max
    ) / (
        1 - (N - 1) * w_max
    ) ** 2 * (
        amp_d
        * tau_d
        * (
            tau_d * (2 + (2 - N) * w_max)
            + tau_syn * (2 - (N - 2) * w_max - (N - 1) * w_max * w_max)
        )
    ) / (
        tau_syn + (1 + w_max) * tau_d
    ) / (
        tau_syn + tau_d * (1 - (N - 1) * w_max)
    )
    return interaction_term


def delta_W_th_truncated(N, rho, w, f0=-0.0013, f10=0.0129762, f11=0.0129762):
    ### f-coefficients computed for default parameters in delta_W_th, see mathematica file
    r = rho / (1 - (N - 1) * w)
    dw = f0 * r * r + 2 * r * f10 * w + r * f11 * w * w * (N - 2)
    return dw


class hawkes_network:
    def __init__(self, **kwargs):

        np.random.seed(kwargs["seed"])
        self.N = kwargs["N"]
        self.w_max = kwargs["w_max"]
        self.rho = kwargs["rho"]
        self.f0 = self.rho * np.ones(self.N)
        self.set_f0()
        self.f0sum = self.f0.sum()
        self.f0prob = self.f0 / self.f0sum
        self.mu = kwargs["mu"]
        self.f1 = np.zeros(self.N)
        ########## symmetric STDP ################
        self.tau_p = kwargs["tau_p"]
        self.tau_d = kwargs["tau_d"]
        self.amp_p = kwargs["amp_p"]
        self.amp_d = kwargs["amp_d"]
        self.traces_spike_update_sym = np.array(
            [self.amp_p, self.amp_p, self.amp_d, self.amp_d]
        )
        self.stdp_traces = np.zeros((self.N, 4))
        ########## for asymmetric STDP ###########
        self.traces_spike_update_asym = np.array([self.amp_p, self.amp_d])
        self.stdp_traces_asym = np.zeros((self.N, 2))
        ########## weight and adjacency matrices ###
        if kwargs["init_W"] == "assemblies":
            self.W = assembly_weight_matrix(
                self.N, kwargs["num_assemblies"], kwargs["assembly_size"], self.w_max
            )
        elif kwargs["init_W"] == "intertwined":
            self.W = intertwined_weight_matrix(
                self.N, self.w_max, kwargs["num_assemblies"], kwargs["swaps"]
            )
        elif kwargs["init_W"] == "random":
            self.W = rand_init_weight_matrix(
                self.N, self.w_max, scale=kwargs["init_scale"]
            )
        self.Adja = np.ones((self.N, self.N))
        np.fill_diagonal(self.Adja, 0)
        ##########################################
        self.t = 0.0
        self.spikes = [[] for _ in range(self.N)]

    def set_f0(self):
        self.f0sum = self.f0.sum()
        self.f0prob = self.f0 / self.f0sum

    def next_spike(self):
        # compute next spiketime: time dt0 until first spontaneous spike
        dt0 = np.random.exponential(1 / self.f0sum)
        f1sum = self.f1.sum()  # time dt1 until first evoked spike
        if f1sum > 100000:
            raise NameError("Activity exploded")
        xunif = np.random.random_sample()
        if xunif >= 1 - np.exp(-f1sum):
            dt1 = dt0 + 1  # dt1 is actually infinity
        else:
            dt1 = CDFm1(f1sum, xunif)  # dt1 is finite
        if dt0 < dt1:  # determine spiking neuron: if next spontaneous spike is first
            dt = dt0
            # distribute it among neurons according to their spont. rates
            spiking_neuron = np.random.choice(self.N, 1, p=self.f0prob)[0]
        else:  # if next evoked spike is first
            dt = dt1
            f1prob = self.f1 / f1sum
            # distribute it among neurons according to their evoked rates
            spiking_neuron = np.random.choice(self.N, 1, p=f1prob)[0]
        self.t += dt
        # neuron frequencies
        self.f1 = self.f1 * np.exp(-dt) + self.W[:, spiking_neuron]

        self.stdp_traces = stdp_prop_sym(dt, self.tau_p, self.tau_d) * self.stdp_traces
        self.stdp_traces[spiking_neuron, :] += self.traces_spike_update_sym

        return dt, spiking_neuron

    def stdp(self, spiking_neuron):
        dW = np.zeros((self.N, self.N))

        dW[spiking_neuron, :] = (
            self.Adja[spiking_neuron, :]
            * (self.stdp_traces[:, 0] + self.stdp_traces[:, 2])
            * self.mu
        )
        dW[:, spiking_neuron] = (
            self.Adja[:, spiking_neuron]
            * (self.stdp_traces[:, 1] + self.stdp_traces[:, 3])
            * self.mu
        )

        return dW

    def weight_update(self, dW, external=0):

        dW_internal = dW[: self.N - external, : self.N - external]
        W_internal = self.W[: self.N - external, : self.N - external]

        W_internal += dW_internal

        # weight restriction
        W_internal[W_internal > self.w_max] = self.w_max
        W_internal[W_internal < 0.0] = 0.0

    def simulate(self, runtime, warmup=False, external=0, tqdmsilence=False):
        tend = self.t + runtime
        #        if not warmup:
        pbar = tqdm(total=runtime, position=0, disable=tqdmsilence or warmup)
        while self.t < tend:
            dt, spiking_neuron = self.next_spike()
            if not warmup:
                dW = self.stdp(spiking_neuron)
                self.weight_update(dW, external)
                pbar.update(dt)
            self.spikes[spiking_neuron].append(self.t)
        pbar.close()

    def track_dw(self, runtime):
        Wtrack = np.zeros((self.N, self.N))
        tend = self.t + runtime
        while self.t < tend:
            _, spiking_neuron = self.next_spike()
            Wtrack += self.stdp(spiking_neuron)
        return Wtrack


def sizes_sim(params):
    network = hawkes_network(**params)
    network.simulate(10000, warmup=True)
    network.simulate(params["tend"], tqdmsilence=True)
    L, W_ = cluster(network.W)
    sizes_louvain = assembly_sizes(L)
    ranges = [int(np.sum(sizes_louvain[:i])) for i in range(len(sizes_louvain) + 1)]
    wsums = [
        np.sum(W_[ranges[i] : ranges[i + 1], ranges[i] : ranges[i + 1]])
        for i in range(len(sizes_louvain))
    ]
    sizes_wsum = [
        0.5 + 0.5 * np.sqrt(1 + 4 * np.sum(wsum / params["w_max"])) for wsum in wsums
    ]
    return (sizes_louvain, sizes_wsum)

