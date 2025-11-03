using Graphs
using SimpleWeightedGraphs
using ProgressBars

"""
    makeShortestPathTree(G::SimpleWeightedGraph, source::Int; target_vertices::Union{Symbol, AbstractVector}=:all)

Construct a shortest path tree from a source vertex in a weighted graph.

# Arguments
- `G::SimpleWeightedGraph`: The input graph with edge weights.
- `source::Int`: The source vertex for the shortest paths.

# Keywords
- `target_vertices::Union{Symbol, AbstractVector}=:all`: List of target vertices. Use `:all` for all vertices.

# Returns
- `SPT::SimpleWeightedGraph`: The shortest path tree.
"""

function unionOfShortestPathTrees(G::SimpleWeightedGraph, anchors, verbose::Bool=false)
    
    # Compute shortest paths (distances and parent vertices) using Dijkstra's algorithm
    # Dijkstra's returns an array of parents. If parent[i] = j, then j is the predecessor of i on a shortest path from source.

	sources = Int[]
	destinations = Int[]
	weights = Real[]

	unique_edges = Set{Tuple{Int, Int}}()

    if verbose
        iter = ProgressBar(anchors)
    else
        iter = anchors
    end
    
    for source in iter
        parents = dijkstra_shortest_paths(G, source).parents

        for target in 1:nv(g)
            if target != source && parents[target] != 0
                u, v = parents[target], target
                if u > v
                    u, v = v, u # Ensure the edge is added in the correct order (smaller vertex first) for canonical representation in the set	
                end

				if !((u,v) in unique_edges)
	                push!( sources, u)
	                push!( destinations, v )
	                push!( weights, g.weights[u,v] )
					push!( unique_edges, (u,v) )
				end
            end
        end
    end
    return SimpleWeightedGraph( sources, destinations, weights )
end

function unionOfAllShortestPathTrees(G::SimpleWeightedGraph, verbose::Bool=false)

    # Determine source vertices
    return unionOfShortestPathTrees(G, 1:nv(G), verbose ) 
end

g = SimpleWeightedGraph(4)
add_edge!(g, 1, 2, 0.5);
add_edge!(g, 1, 3, 0.8);
add_edge!(g, 2, 3, 2);
add_edge!(g, 3, 4, 0.2); 
add_edge!(g, 1, 4, 0.7);

spt = unionOfShortestPathTrees(g, [1, 2, 3], true)

mb = unionOfAllShortestPathTrees(g, true)

for e in edges(spt)
    println("Edge: $(e.src) -> $(e.dst) with weight $(e.weight)")
end