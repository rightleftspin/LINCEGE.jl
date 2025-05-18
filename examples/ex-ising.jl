using JLD
using NLCE

order = 8
num_sites, bond_lists, multiplicities = load("./outputs/Kagome_Lattice_$(order).jld", "num_sites"),
load("./outputs/Kagome_Lattice_$(order).jld", "bond_lists"),
load("./outputs/Kagome_Lattice_$(order).jld", "multiplicities")

B = 0
couplings =  [1, 5e-2]
temperature = range(0, 10, length=100)

obs = NLCE.ising_observables(num_sites, bond_lists, multiplicities, temperature, B, couplings, 3, 6)

save("./outputs/kagome_lattice_obs_$(order).jld", "temp", temperature, "obs", obs)
