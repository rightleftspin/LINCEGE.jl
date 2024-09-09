What I need to do

Create Lattice, graph object
    Lattice Object, currently only support for infinite lattice

Lattice Generators creates graph object which is not a predefined thing, its just an
output from the generator function

break lattice

Functions:
    Generate Lattice: Directions -> Graph Object
    -- Tile Graph: Graph Object, Finite or Infinite -> Array of graph objects
    Count Graphs: Array of Graphs -> dict of graphs with counts
    Count Subgraphs: Array of Graphs -> dict of graphs with counts
    NLCE Summation: Dict of graphs with counts -> dict of graphs with counts

Features:
    Simple NLCE using bravais lattice
    Simple NLCE using unit cell
    Tiling NLCE using cell
    Tiling NLCE using unit cell and cell


Pipeline:
    
    Step 1: Break Lattice
        Input: Lattice
        Output: Array of Lattice Broken Down
    Step 2: Filter Array of Clusters
        Input: Array of Clusters
        Output: 
