### Global parameters ###

seconds_per_unit_time = 0.01

#########################

pars_spont = {
    "tau_p": 2.5,
    "tau_d": 5.0,
    "amp_p": 0.08,
    "amp_d": -0.0533,
    "rho": 0.0015,
    "N": 50,
    "w_max": 0.05,
    "mu": 0.07,
    "seed": None,
    "tend": 50_000_000,
    "r_in": 0.04,
    "w_in": 0.05,
    "init_W": "random",
    "init_scale": 0.2,
}

pars_avg_dw = {
    "tau_p": 2.5,
    "tau_d": 5.0,
    "amp_p": 0.08,
    "amp_d": -0.0533,
    "rho": 0.0015,
    "N": 50,
    "w_max": 0.05,
    "mu": 0.07,
    "seed": None,
    "tend": 50_000_000,
    "init_W": None,
}

pars_learn = {
    "tau_p": 3.5,
    "tau_d": 5.0,
    "amp_p": 0.08,
    "amp_d": -0.065,
    "rho": 0.0015,
    "rho_ext": 0.0418,
    "N": 81,
    "w_max": 0.026,
    "w_ext": 0.26,
    "mu": 0.07,
    "seed": None,
    "assembly_size": 20,
    "inputs": 1,
    "t_ON": 18_000,
    "t_OFF": 10_000_000,
    "init_W": "random",
    "init_scale": 0.1,
}


pars_drift = {
    "tau_p": 2.5,
    "tau_d": 5.0,
    "amp_p": 0.08,
    "amp_d": -0.0533,
    "rho": 0.002,
    "N": 72,
    "w_max": 0.056,
    "mu": 0.148,
    "seed": None,
    "T1": 50_000_000,
    "T2": 50_000_000,
    "init_W": "random",
    "init_scale": 0.25,
}


pars_drift2 = {
    "tau_p": 2.5,
    "tau_d": 5.0,
    "amp_p": 0.08,
    "amp_d": -0.0533,
    "rho": 0.0015,
    "rho_small": 0.0003,
    "N": 120,
    "w_max": 0.024,
    "mu": 0.05,
    "seed": None,
    "t_switch": 30_000_000,
    "p_switch": 0.03,
    "init_W": "assemblies",
    "num_assemblies": 6,
    "assembly_size": 20,
}

pars_sizes = {
    "tau_p": 2.5,
    "tau_d": 5.0,
    "amp_p": 0.08,
    "amp_d": -0.0533,
    "rho": 0.0015,
    "N": 150,
    "mu": 0.04,
    "seed": None,
    "tend": 150_000_000,
    "init_W": "random",
    "init_scale": 0.2,
}


pars_intertwined = {
    "seconds_per_unit_time": 0.01,
    "tau_p": 2.6,
    "tau_d": 6.5,
    "amp_p": 0.08,
    "amp_d": -0.042,
    "rho": 0.0015,
    "w_max": 0.018,
    "N": 190,
    "num_assemblies": 20,
    "swaps": 0,
    "mu": 0.017,
    "seed": None,
    "t_eq": 20_000_000,
    "n_sims": 900,
    "t_sim": 100_000,
    "init_W": "intertwined",
}

pars_avg_dw = {
    "tau_p": 2.5,
    "tau_d": 5.0,
    "amp_p": 0.08,
    "amp_d": -0.0533,
    "rho": 0.0015,
    "N": 50,
    "w_max": 0.05,
    "mu": 0.07,
    "seed": None,
    "tend": 50_000_000,
    "init_W": None,
}

pars_overlap = {
    "tau_p": 2.5,
    "tau_d": 5.0,
    "amp_p": 0.08,
    "amp_d": -0.0533,
    "rho": 0.0015,
    "rho_small": 0.0001,
    "N": 60,
    "w_max": 0.024,
    "mu": 0.045,
    "seed": None,
    "t_end": 100_000_000,
    "init_W": "assemblies",
    "num_assemblies": 3,
    "assembly_size": 20,
}


pars_sparse = {
    "tau_p": 2.5,
    "tau_d": 5.0,
    "amp_p": 0.08,
    "amp_d": -0.0533,
    "rho": 0.0015,
    "N": 50,
    "w_max": 0.05,
    "mu": 0.07,
    "seed": None,
    "tend": 20_000_000,
    "init_W": None,
    "density": 0.8,
}

pars_input_strength = {
    "tau_p": 3.5,
    "tau_d": 5.0,
    "amp_p": 0.08,
    "amp_d": -0.066,
    "rho": 0.0015,
    "N": 50,
    "N_target": 20,
    "w_max": 0.026,
    "mu": 0.01,
    "seed": None,
    "r_in": 0.04,
    "w_in": 0.05,
    "init_W": None,
}
