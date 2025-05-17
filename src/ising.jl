"""
Utility functions for testing the Ising model
"""
function ising_model(num_sites, bonds, B, couplings)
    J, mu = couplings
    num_spins = 2 ^ num_sites
    energies = zeros(num_spins)
    magnetizations = zeros(num_spins)

    for spin_config in 0:(num_spins - 1)
        spins = (2 * digits(spin_config, base=2, pad=num_sites)) .- 1
        energies[spin_config + 1] -= B * mu * sum(spins)
        for bond in bonds
          if spins[bond[1]] == spins[bond[2]]
            energies[spin_config + 1] += J
          else
            energies[spin_config + 1] -= J
          end
        end
        magnetizations[spin_config + 1] = sum(spins)
    end

    energies, magnetizations
end

function ising_observables(sites_per_cluster, bond_lists, multiplicities, temperatures, B, couplings)
  resums = 0
  # Energy, Entropy, Specific Heat, Magnetization
  # Extra resum multiplicities for euler and wynn
  # Properties: property, order, temperature
  properties = zeros(4, length(temperature), length(multiplicities[1]) + resums)

  for (num_sites, bond_list, mult) in zip(sites_per_cluster, bond_lists, multiplicities)
    eigs, mags = ising_model(num_sites, bond_list, B, couplings)

    exp_energy_temp_matrix = exp.(-transpose(eigs) ./ temperatures)
    partition_function = sum(exp_energy_temp_matrix, dims=2)
    avg_energy = (exp_energy_temp_matrix * eigs) ./ partition_function

    mult_plus_resum = transpose(append!(mult, repeat([0], resums)))
    properties[1, :, :] += mult_plus_resum .* avg_energy
    properties[2, :, :] += mult_plus_resum .* (log.(partition_function) + (avg_energy ./ temperatures))
    properties[3, :, :] += mult_plus_resum .* ((((exp_energy_temp_matrix * (eigs .^ 2)) ./ partition_function) - (avg_energy .^ 2)) ./ (temperatures .^ 2))
    properties[4, :, :] += mult_plus_resum .* ((exp_energy_temp_matrix * mags) ./ partition_function)

  end
      properties
end

using Plots
using NLCE

square_lattice = Dict("Basis" => [[0, 0]], "Primitive Vectors" => [[1, 0], [0, 1]])

neighbors = [1]
max_order = 8

square_nlce_bundle = NLCE.SiteExpansionBundle(
    square_lattice["Basis"],
    square_lattice["Primitive Vectors"],
    neighbors,
    max_order,
    NLCE.isomorphic_pruning,
)

square_lattice_cluster_info = NLCE.lattice_constants!(
    square_nlce_bundle,
    length(NLCE.start(square_nlce_bundle)),
)

NLCE.subclusters!(square_nlce_bundle, false)

final_weights = NLCE.final_clusters(square_nlce_bundle)

num_sites = []
bond_lists = []
multiplicities = []
for (cluster, mults) in final_weights
    push!(num_sites, NLCE.nv(cluster))
    push!(bond_lists, NLCE.weighted_edge_list(cluster))
    push!(multiplicities, mults)
end


B = 4
couplings = [0.5, 5e-2]
temperature = range(0, 10, length=100)

obs = ising_observables(num_sites, bond_lists, multiplicities, temperature, B, couplings)

plot(temperature, obs[1, :, end-3:end])
savefig("avg_energy_ising.pdf")

plot(temperature, obs[2, :, end-3:end])
savefig("entropy_ising.pdf")

plot(temperature, obs[3, :, end-3:end])
savefig("cv_ising.pdf")

plot(temperature, obs[4, :, end-3:end], ylimits=(0,1), xlimits=(0,5))
savefig("magnetization_ising.pdf")
