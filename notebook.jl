### A Pluto.jl notebook ###
# v0.20.8

using Markdown
using InteractiveUtils

# ╔═╡ 2ce47dbc-dac3-4613-8873-e1a347d0397a
using Graphs, SimpleWeightedGraphs

# ╔═╡ d0101353-a6c0-4b5d-b732-a321b0f7552e
using DataStructures # For `Set` which is similar to Python's set

# ╔═╡ 7a90d16a-6c6d-4dde-838e-2aa7527897c0
using ProgressMeter # For `@showprogress` similar to `tqdm`

# ╔═╡ 83ec0c60-1b3f-4e27-a746-0d6247e88212
using ProgressBars

# ╔═╡ 16eb8e1c-0f4f-4513-bf7d-26c6a3bbf06a
include("metric_sparsification.jl") 

# ╔═╡ 351282fe-742b-11f0-29a6-0b7b1d1735a7
# ╠═╡ disabled = true
#=╠═╡
import Graphs
  ╠═╡ =#

# ╔═╡ 3e41c919-3dc2-4c04-b612-21a7571709a6
begin
	g = SimpleWeightedGraph(4)
	add_edge!(g, 1, 2, 0.5);
	add_edge!(g, 1, 3, 0.8);
	add_edge!(g, 2, 3, 2);
	add_edge!(g, 3, 4, 0.2); 
	add_edge!(g, 1, 4, 0.7);
end

# ╔═╡ 4fe58a45-8147-4803-b111-0a748f08e80d
for e in edges(g)
	print(e)
end

# ╔═╡ dfafb223-7d01-4d8e-a15b-1c5eb8be941f
dist = dijkstra_shortest_paths(g, 1)

# ╔═╡ 703488b0-287c-442f-9b04-2dd708c1a546
g.weights[1,2]

# ╔═╡ 7a3ade7b-e612-4f31-abfa-8e3479cd2a33
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
function makeShortestPathTree(G::SimpleWeightedGraph, source::Int; target_vertices::Union{Symbol, AbstractVector}=:all)
    
    # Determine target vertices
    if target_vertices == :all
        target_vertices_list = 1:nv(G)
    else
        target_vertices_list = target_vertices
    end

    # Compute shortest paths (distances and parent vertices) using Dijkstra's algorithm
    # Dijkstra's returns an array of parents. If parent[i] = j, then j is the predecessor of i on a shortest path from source.
    dists = dijkstra_shortest_paths(G, source).dists
    parents = dijkstra_shortest_paths(G, source).parents

    SPT = SimpleWeightedGraph(nv(G))
    
    # Collect unique edges that form the shortest paths
    unique_edges = Set{Tuple{Int, Int}}()

    for target in target_vertices_list
        current = target
        while current != source && parents[current] != 0
            # Ensure the edge is added in the correct order (smaller vertex first) for canonical representation in the set
            u, v = parents[current], current
            if u > v
                u, v = v, u
            end
            push!(unique_edges, (u, v))
            current = parents[current]
        end
    end

    # Add edges and their weights to the Shortest Path Tree graph
    for (u, v) in unique_edges
        if has_edge(G, u, v)
            add_edge!(SPT, u, v, get_weight(G, u, v))
        end
    end

    return SPT
end


# Helper function to add all edges from one graph to another

# ╔═╡ 6bd3e16d-d3b8-4cd8-b2ba-c0880763d097
# Helper function to add all edges from one graph to another
function add_graph_edges!(G_dest::SimpleWeightedGraph, G_source::SimpleWeightedGraph)
    for edge in edges(G_source)
        u, v = src(edge), dst(edge)
        weight = get_weight(G_source, u, v)
        if !has_edge(G_dest, u, v)
            add_edge!(G_dest, u, v, weight)
        end
    end
    return G_dest
end


# ╔═╡ 391b7efc-65fb-4138-a253-54d79319c6cb
"""
    union_shortest_path_trees(G::SimpleWeightedGraph, source_vertices::Union{Symbol, AbstractVector}=:all, target_vertices::Union{Symbol, AbstractVector}=:all; verbose::Bool=false)

Compute the union of shortest path trees from source vertices to target vertices.

# Arguments
- `G::SimpleWeightedGraph`: A graph with edges weighted as distances.
- `source_vertices::Union{Symbol, AbstractVector}`: List of source vertices. Use `:all` for all vertices.
- `target_vertices::Union{Symbol, AbstractVector}`: List of target vertices. Use `:all` for all vertices.

# Keywords
- `verbose::Bool=false`: Whether to display a progress bar.

# Returns
- `G_union::SimpleWeightedGraph`: The union of shortest path trees.
"""
function union_shortest_path_trees(G::SimpleWeightedGraph, source_vertices::Union{Symbol, AbstractVector}=:all, target_vertices::Union{Symbol, AbstractVector}=:all; verbose::Bool=false)
    
    # Initialize an empty graph with the same number of vertices as G
    G_union_spt = SimpleWeightedGraph(nv(G))

    # Determine source vertices
    if source_vertices == :all
        source_vertices_iter = 1:nv(G)
    else
        source_vertices_iter = source_vertices
    end

    # Iterate through source vertices with optional progress bar
    if verbose
        p = Progress(length(source_vertices_iter), 1, "Computing SPTs and Union...")
        for source in source_vertices_iter
            SPT_source = makeShortestPathTree(G, source, target_vertices=target_vertices)
            G_union_spt = add_graph_edges!(G_union_spt, SPT_source)
            next!(p)
        end
    else
        for source in source_vertices_iter
            SPT_source = makeShortestPathTree(G, source, target_vertices=target_vertices)
            G_union_spt = add_graph_edges!(G_union_spt, SPT_source)
        end
    end
    
    return G_union_spt
end

# ╔═╡ 496cdf95-608f-4423-8143-ea8df054889b
spt_1 = makeShortestPathTree( g, 1 )

# ╔═╡ fa998c17-ba7c-44e8-9e5c-163e02a1d023
spt_2 = makeShortestPathTree( g, 2 )

# ╔═╡ 78f83336-0e75-4ae0-95e5-60daf1a18673
begin
	source = 2 
	parents = dijkstra_shortest_paths(g, source).parents
end

# ╔═╡ c9acb016-c0eb-477b-ba24-c1521ada3f51
unique_edges = Set{Tuple{Int, Int, Real}}()

# ╔═╡ 841541f9-f376-47f5-9465-97a9a374bbc3
begin
	sources = Int[];
	destinations = Int[];
	weights = Real[]
end

# ╔═╡ d5df8d07-7908-4e13-b12e-35e749c696df
for target in 1:nv(g)
    if target != source && parents[target] != 0
        u, v = parents[target], target
        if u > v
            u, v = v, u # Ensure the edge is added in the correct order (smaller vertex first) for canonical representation in the set	
    	end
		push!( sources, u)
		push!( destinations, v )
		push!( weights, g.weights[u,v] )

        push!(unique_edges, (u, v, g.weights[u,v]))
    end
end

# ╔═╡ c29e4fc4-ddc2-49fb-ae6a-abd7a9ba23ad
sources

# ╔═╡ afa86854-481e-4ce9-bb82-b679d4d1f251
destinations

# ╔═╡ 6839920f-8f11-4afa-b874-179f9806ba5b
weights

# ╔═╡ dcfc0862-bdbd-48d8-a202-ab8f116290fd
SPT = SimpleWeightedGraph( sources, destinations, weights )

# ╔═╡ 344ffa0e-5e12-45a8-8e7d-89fc410c40d0
for e in edges(SPT)
	print(e)
end

# ╔═╡ fc23fd5c-69ea-4604-964d-d147ff0b7890
G_mb = union_shortest_path_trees( g )

# ╔═╡ 59773f89-f761-46c0-a10e-31a38692c8de
for e in edges(G_mb)
	print(e)
end

# ╔═╡ 3f3a0255-800d-4013-a9dd-8cd9750229c0
1:5

# ╔═╡ fe9a09a3-25e3-4391-9ffd-ee5b8f9b3e12
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


# ╔═╡ 98d76ea0-bfb5-4076-b778-0c8ca139bbc6
essai = unionOfShortestPathTrees(g, 1:4, true)

# ╔═╡ d0328763-ea75-4e73-bea2-bf42f37655d6
for e in edges(essai)
	print(e)
end

# ╔═╡ 6be52c38-339c-4e5c-89a1-ca2194567828
for e in edges(g)
	print(e)
end

# ╔═╡ 0dd88a7e-50b5-4c38-bf5b-1ed516484eac
for i in ProgressBar(1:2) #wrap any iterator
	a = 1
       end

# ╔═╡ b36042fd-5299-4a85-8de5-72cbd1bbca97


# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DataStructures = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
Graphs = "86223c79-3864-5bf0-83f7-82e725a168b6"
ProgressBars = "49802e3a-d2f1-5c88-81d8-b72133a6f568"
ProgressMeter = "92933f4c-e287-5a05-a399-4b506db050ca"
SimpleWeightedGraphs = "47aef6b3-ad0c-573a-a1e2-d07658019622"

[compat]
DataStructures = "~0.18.22"
Graphs = "~1.13.0"
ProgressBars = "~1.5.1"
ProgressMeter = "~1.10.4"
SimpleWeightedGraphs = "~1.5.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.6"
manifest_format = "2.0"
project_hash = "280fa2038d8c3914760117e364b3ad868330ad3f"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "d57bd3762d308bded22c3b82d033bff85f6195c6"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.4.0"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "8ae8d32e09f0dcf42a36b90d4e17f5dd2e4c4215"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.16.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "4e1fe97fdaed23e9dc21d4d664bea76b65fc50a0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.22"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
version = "1.11.0"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "c5abfa0ae0aaee162a3fbb053c13ecda39be545b"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.13.0"

[[deps.Inflate]]
git-tree-sha1 = "d1b1b796e47d94588b3757fe84fbf65a5ec4a80d"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.5"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OrderedCollections]]
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.ProgressBars]]
deps = ["Printf"]
git-tree-sha1 = "b437cdb0385ed38312d91d9c00c20f3798b30256"
uuid = "49802e3a-d2f1-5c88-81d8-b72133a6f568"
version = "1.5.1"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "13c5103482a8ed1536a54c08d0e742ae3dca2d42"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.10.4"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"
version = "1.11.0"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "be8eeac05ec97d379347584fa9fe2f5f76795bcb"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.5"

[[deps.SimpleWeightedGraphs]]
deps = ["Graphs", "LinearAlgebra", "Markdown", "SparseArrays"]
git-tree-sha1 = "3e5f165e58b18204aed03158664c4982d691f454"
uuid = "47aef6b3-ad0c-573a-a1e2-d07658019622"
version = "1.5.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.11.0"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "cbea8a6bd7bed51b1619658dec70035e07b8502f"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.14"

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

    [deps.StaticArrays.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.7.0+0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"
"""

# ╔═╡ Cell order:
# ╠═351282fe-742b-11f0-29a6-0b7b1d1735a7
# ╠═2ce47dbc-dac3-4613-8873-e1a347d0397a
# ╠═3e41c919-3dc2-4c04-b612-21a7571709a6
# ╠═4fe58a45-8147-4803-b111-0a748f08e80d
# ╠═dfafb223-7d01-4d8e-a15b-1c5eb8be941f
# ╠═703488b0-287c-442f-9b04-2dd708c1a546
# ╠═d0101353-a6c0-4b5d-b732-a321b0f7552e
# ╠═7a90d16a-6c6d-4dde-838e-2aa7527897c0
# ╠═83ec0c60-1b3f-4e27-a746-0d6247e88212
# ╠═391b7efc-65fb-4138-a253-54d79319c6cb
# ╠═7a3ade7b-e612-4f31-abfa-8e3479cd2a33
# ╠═6bd3e16d-d3b8-4cd8-b2ba-c0880763d097
# ╠═496cdf95-608f-4423-8143-ea8df054889b
# ╠═fa998c17-ba7c-44e8-9e5c-163e02a1d023
# ╠═78f83336-0e75-4ae0-95e5-60daf1a18673
# ╠═c9acb016-c0eb-477b-ba24-c1521ada3f51
# ╠═841541f9-f376-47f5-9465-97a9a374bbc3
# ╠═d5df8d07-7908-4e13-b12e-35e749c696df
# ╠═c29e4fc4-ddc2-49fb-ae6a-abd7a9ba23ad
# ╠═afa86854-481e-4ce9-bb82-b679d4d1f251
# ╠═6839920f-8f11-4afa-b874-179f9806ba5b
# ╠═dcfc0862-bdbd-48d8-a202-ab8f116290fd
# ╠═344ffa0e-5e12-45a8-8e7d-89fc410c40d0
# ╠═fc23fd5c-69ea-4604-964d-d147ff0b7890
# ╠═59773f89-f761-46c0-a10e-31a38692c8de
# ╠═3f3a0255-800d-4013-a9dd-8cd9750229c0
# ╠═fe9a09a3-25e3-4391-9ffd-ee5b8f9b3e12
# ╠═98d76ea0-bfb5-4076-b778-0c8ca139bbc6
# ╠═d0328763-ea75-4e73-bea2-bf42f37655d6
# ╠═6be52c38-339c-4e5c-89a1-ca2194567828
# ╠═0dd88a7e-50b5-4c38-bf5b-1ed516484eac
# ╠═16eb8e1c-0f4f-4513-bf7d-26c6a3bbf06a
# ╠═b36042fd-5299-4a85-8de5-72cbd1bbca97
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
