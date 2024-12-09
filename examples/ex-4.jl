"""
Example generating the data for a kagome lattice and 
writing to a file meant to be read by fortran code
"""
module ex1

using NLCE

# Set the basis, since there is only one atom, it is at [0, 0]
basis = [[0, 0, 0], [0, 1 / 4, 1 / 4], [1 / 4, 0, 1 / 4], [1 / 4, 1 / 4, 0]]

# Choose the primitive vectors, there are two on a square lattice
primitive_vec = [[0, 1 / 2, 1 / 2], [ 1/ 2, 0, 1/2], [1 / 2, 1 / 2, 0]]

# Choosing nearest neighbors (within distance 1 from each other)
neighborhood = [1 / 2]

# Setting the maximum order
max_order = 4

# Generating all the clusters using this information
nlce_clusters = simple_NLCE(basis, primitive_vec, neighborhood, max_order)

# Writing all the files to the corresponding folder, creating the folder
# if it does not exist
filepath = "examples/outputs/ex-4/pyrochlore_nn"
mkpath(filepath)
filename = filepath * "/pyrochlore_nn"

# Write all the files in the "fortran" format
write_to_file_fortran(nlce_clusters, filename, max_order)

end
