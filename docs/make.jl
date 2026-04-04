using Lincege
using Documenter

DocMeta.setdocmeta!(Lincege, :DocTestSetup, :(using Lincege); recursive=true)

makedocs(;
        modules=[Lincege],
        authors="Pranav Seetharaman <pranav@myrdd.info> and contributors",
        sitename="Lincege.jl",
        format=Documenter.HTML(;
                canonical="https://rightleftspin.github.io/Lincege.jl",
                edit_link="main",
                assets=String[],
        ),
        pages=[
                "Home" => "index.md",
                "Examples" => [
                        "Square Lattice" => "examples/square_lattice.md",
                        "Kagome Lattice" => "examples/kagome_lattice.md",
                        "Pyrochlore Lattice Unit Cell" => "examples/pyrochlore_unit_cell_expansion.md",
                        "Square Lattice Cluster" => "examples/square_cluster.md",
                ],
        ],
)

deploydocs(;
        repo="github.com/rightleftspin/Lincege.jl",
        devbranch="main",
)
