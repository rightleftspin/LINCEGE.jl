
using NLCE, Test

using Base.Iterators

@testset verbose = true "NLCE" begin

    # A few different lattices to test on
    square_lattice = Dict("Basis" => [[0, 0]], "Primitive Vectors" => [[1, 0], [0, 1]])

    triangular_lattice =
        Dict("Basis" => [[0, 0]], "Primitive Vectors" => [[0.5, sqrt(3) / 2], [1, 0]])

    kagome_lattice = Dict(
        "Basis" => [[0, 0], [1, 0], [1 / 2, sqrt(3) / 2]],
        "Primitive Vectors" => [[2, 0], [1, sqrt(3)]],
    )

    @testset "Bundles" begin
        # Here I will test a variety of lattices at nearest neighbor up till order 6
        # to see if the Cluster struct works. Most issues usually present themselves
        # by order 5.

        # Counts according to Rigol, Bryant and Singh 2007
        num_topo_clusters_square = [1, 1, 1, 3, 4, 10, 19, 51, 112, 300]
        sum_lc_clusters_square = [1, 2, 6, 19, 63, 216, 760, 2725, 9910, 36446]

        num_topo_clusters_triangular = [1, 1, 2, 4, 8, 22, 54, 156, 457, 1424]
        sum_lc_clusters_triangular = [1, 3, 11, 44, 186, 814, 3652, 16689, 77359, 362671]

        num_topo_clusters_kagome = [1, 1, 2, 2, 4, 7, 12, 22, 45, 88]
        sum_lc_clusters_kagome = [1, 2, 14//3, 12, 33, 281//3, 272, 805, 2420, 7358]

        neighbors = [1]
        max_order = 4

        square_nlce_bundle = NLCE.SiteExpansionBundle(
            square_lattice["Basis"],
            square_lattice["Primitive Vectors"],
            neighbors,
            max_order,
            NLCE.isomorphic_pruning
        )

        triangular_nlce_bundle = NLCE.SiteExpansionBundle(
            triangular_lattice["Basis"],
            triangular_lattice["Primitive Vectors"],
            neighbors,
            max_order,
            NLCE.isomorphic_pruning
        )

        kagome_nlce_bundle = NLCE.SiteExpansionBundle(
            kagome_lattice["Basis"],
            kagome_lattice["Primitive Vectors"],
            neighbors,
            max_order,
            NLCE.isomorphic_pruning
        )

        square_lattice_cluster_info = NLCE.lattice_constants!(square_nlce_bundle)
        triangular_lattice_cluster_info = NLCE.lattice_constants!(triangular_nlce_bundle)
        kagome_lattice_cluster_info = NLCE.lattice_constants!(kagome_nlce_bundle)

        @test sum([mult for (hash, (_, mult, _, _, _)) in square_lattice_cluster_info]) == sum(sum_lc_clusters_square[1:max_order])
        @test length(square_lattice_cluster_info) == sum(num_topo_clusters_square[1:max_order])

        @test sum([mult for (hash, (_, mult, _, _, _)) in triangular_lattice_cluster_info]) == sum(sum_lc_clusters_triangular[1:max_order])
        @test length(triangular_lattice_cluster_info) == sum(num_topo_clusters_triangular[1:max_order])

        @test sum([mult for (hash, (_, mult, _, _, _)) in kagome_lattice_cluster_info]) == sum(sum_lc_clusters_kagome[1:max_order])
        @test length(kagome_lattice_cluster_info) == sum(num_topo_clusters_kagome[1:max_order])

        # TODO: Write tests for subcluster multiplicities

    end

    @testset "Util" begin
        # Testing various utility functions
        test_adj_list = [[2, 3], [1, 4], [1, 4], [2, 3]]
        test_connections = [[1, 2], [2, 3], [4, 5], [5, 6]]

        test_adj_matrices = zeros(Int, 3, 6, 6)
        test_adj_matrices[1, :, :] = [1 1 0 0 0 0;
                                      1 1 1 0 0 0;
                                      0 1 1 0 0 0;
                                      0 0 0 1 1 0;
                                      0 0 0 1 1 1;
                                      0 0 0 0 1 1;]

        test_adj_matrices[2, :, :] = [1 1 0 0 0 0;
                                      3 1 1 0 0 0;
                                      0 3 1 0 0 0;
                                      0 0 0 1 1 0;
                                      0 0 0 3 1 1;
                                      0 0 0 0 3 1;]

        test_adj_matrices[3, :, :] = [0 1 0 0 0 0;
                                      1 0 2 0 0 0;
                                      0 2 0 0 0 0;
                                      0 0 0 0 3 0;
                                      0 0 0 3 0 4;
                                      0 0 0 0 4 0;]

        ans1 = zeros(Int, 3, 2, 2)
        ans1[1, :, :] = [1 1;
                         1 1;]
        ans1[2, :, :] = [1 1;
                         3 1;]
        ans1[3, :, :] = [0 1;
                         1 0;]

        @test NLCE.reindex_adj_list(test_adj_list, [4]) == [[]]
        @test NLCE.reindex_adj_list(test_adj_list, [1, 2, 3]) == [[2, 3], [1], [1]]
        @test NLCE.reindex_adj_list(test_adj_list, [1, 2, 4]) == [[2], [1, 3], [2]]

        @test NLCE.reindex_connections(test_connections, [1], [1, 2]) == [[1, 2]]
        @test NLCE.reindex_connections(test_connections, [1, 2, 4], [1, 2, 3, 5, 6]) == [[1, 2], [2, 3], [4, 5]]

        @test NLCE.reindex_adjacency_matrices(test_adj_matrices, [1], [1, 2]) == ans1

    end
end
