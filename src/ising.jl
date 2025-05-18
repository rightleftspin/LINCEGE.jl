"""
Utility functions for testing the Ising model
"""
function ising_model(num_sites, bonds, B, couplings)
    J, mu = couplings
    num_spins = 2 ^ num_sites
    energies = zeros(num_spins)
    magnetizations = zeros(num_spins)

    for spin_config = 0:(num_spins-1)
        spins = (2 * digits(spin_config, base = 2, pad = num_sites)) .- 1
        energies[spin_config+1] -= B * mu * sum(spins)
        for bond in bonds
            if spins[bond[1]] == spins[bond[2]]
                energies[spin_config+1] += J
            else
                energies[spin_config+1] -= J
            end
        end
        magnetizations[spin_config+1] = sum(spins)
    end

    energies, magnetizations
end

function ising_observables(
    sites_per_cluster,
    bond_lists,
    multiplicities,
    temperatures,
    B,
    couplings,
    wynn_cycles,
    euler_start,
)
    resums = 2
    # Energy, Entropy, Specific Heat, Magnetization
    # Extra resum multiplicities for euler and wynn
    # Properties: property, order, temperature
    properties = zeros(4, length(temperatures), length(multiplicities[1]) + resums)

    for (num_sites, bond_list, mult) in zip(sites_per_cluster, bond_lists, multiplicities)
        eigs, mags = ising_model(num_sites, bond_list, B, couplings)

        exp_energy_temp_matrix = exp.(-transpose(eigs) ./ temperatures)
        partition_function = sum(exp_energy_temp_matrix, dims = 2)
        avg_energy = (exp_energy_temp_matrix * eigs) ./ partition_function

        mult_plus_resum = transpose(append!(mult, repeat([0], resums)))
        properties[1, :, :] += mult_plus_resum .* avg_energy
        properties[2, :, :] +=
            mult_plus_resum .* (log.(partition_function) + (avg_energy ./ temperatures))
        properties[3, :, :] +=
            mult_plus_resum .* (
                (
                    ((exp_energy_temp_matrix * (eigs .^ 2)) ./ partition_function) -
                    (avg_energy .^ 2)
                ) ./ (temperatures .^ 2)
            )
        properties[4, :, :] +=
            mult_plus_resum .* ((exp_energy_temp_matrix * mags) ./ partition_function)

    end

    properties[:, :, end - 1] = euler_resummation(properties, euler_start)
    properties[:, :, end] = wynn_resummation(properties, wynn_cycles)
    properties
end
