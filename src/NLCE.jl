module NLCE

include("generation.jl")
include("tiling.jl")
include("tagging.jl")
include("helpers.jl")
include("summation.jl")

include("ising.jl")


# Below is an example of the NLCE process for a square lattice with only nearest neighbors
# up till order 4 it uses a wrapper function for speed, check the helper functions for 
# more info

using AlgebraicNumbers: AlgebraicNumber as AN
using Plots; gr()

## Tasks before next meeting
# do square, triangle and kagome lattice ising models

## Setting up the lattice geometry
basis = [[AN(0), AN(0)]]
primitive_vectors = [[AN(1), AN(0)], [AN(0), AN(1)]]
final_order = 8

xmax = 2
total_temp = collect(0.05:0.05:xmax)
all_energies = []

for order in 7:final_order
    total_energy = zeros(length(total_temp))
    nlce_weights = simple_nlce(basis, primitive_vectors, order)
    for (cluster, mult) in nlce_weights
        ising_energy = ising_energies(cluster)
        total_energy += mult .* energy_solver(total_temp, ising_energy) 
    end
    push!(all_energies, total_energy)
end

plot(total_temp, all_energies, label = collect(7:final_order))
xlims!(0, xmax)
ylims!(-1, 0)
title!("Ising Energy for a Square Lattice")
xlabel!("Temperature")
ylabel!("Energy")
savefig("ising_energy.png")

end
