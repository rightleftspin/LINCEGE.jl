@testset verbose = true "UnitCells" begin

        @testset "Site Expansion Structure" begin
                @test basis_size(square_uc) == 1
                @test dimension(square_uc) == 2
                @test shift_unit_cell(square_uc, [1, 0, 1]) == [1.0, 0.0]
                @test shift_unit_cell(square_uc, [2, -2, 1]) == [2.0, -2.0]
                @test shift_unit_cell(square_uc, [1 2 3; 0 -2 1; 1 1 1]) == [1.0 2.0 3.0; 0.0 -2.0 1.0]
        end

        @testset "Strong Cluster Expansion Structure" begin
                @test basis_size(pyro_exp_uc_uc) == [4]
                @test dimension(pyro_exp_uc_uc) == 3
                @test shift_unit_cell(pyro_exp_uc_uc, [1, -1, 1, 1, 1]) == [1 / 2, 4.5, 1 / 2]
                @test shift_unit_cell(pyro_exp_uc_uc, [1 2 3; -1 1 4; 1 0 -2; 1 1 1; 1 2 4]) == [0.5 6.5 13.5; 4.5 3.5 1.5; 0.5 1.5 4.5]
                @test pyro_exp_uc_uc.translation_labels == [[1, 2, 3, 4]]

        end

        @testset "Weak Cluster Expansion Structure" begin
                @test basis_size(square_cluster_uc) == [4]
                @test dimension(square_cluster_uc) == 2
                @test shift_unit_cell(square_cluster_uc, [1, -1, 1, 1]) == [-0.5, 1.5]
                @test shift_unit_cell(square_cluster_uc, [1 2 3; -1 1 4; 1 1 1; 1 2 4]) == [-0.5 2.5 7.5; 1.5 1.5 -0.5]
                @test square_cluster_uc.translation_labels == [[1, 2, 2, 1]]

        end

        @testset "Visualization" begin
                @test image_unit_cell(square_uc) isa Plots.Plot{Plots.GRBackend}
                @test image_unit_cell(kagome_uc) isa Plots.Plot{Plots.GRBackend}
                @test image_unit_cell(one_fifth_uc) isa Plots.Plot{Plots.GRBackend}
                @test image_unit_cell(ss_uc) isa Plots.Plot{Plots.GRBackend}
                @test image_unit_cell(cube_uc) isa Plots.Plot{Plots.GRBackend}
        end

end # UnitCells
