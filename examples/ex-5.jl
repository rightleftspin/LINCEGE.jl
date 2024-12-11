module ex4

using NLCE

basis = [[0, 0, 0], [1 / 4, 1 / 4, 1 / 4]]

primitive_vec = [[0, 1 / 2, 1 / 2], [1 / 2, 0, 1 / 2], [1 / 2, 1 / 2, 0]]

neighborhood = [sqrt(3) / 4]

# Setting the maximum order
max_order = 8

# Generating all the clusters using this information
nlce_clusters = simple_NLCE(basis, primitive_vec, neighborhood, max_order)

# Writing all the files to the corresponding folder, creating the folder
# if it does not exist
filepath = "examples/outputs/ex-5/diamond_nn"
mkpath(filepath)
filename = filepath * "/diamond_nn"

# Write all the files in the default format
write_to_file(nlce_clusters, filename)

end
