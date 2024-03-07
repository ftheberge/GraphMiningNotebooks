using CGE

edges, weights, vweights, comm, clusters, embed, verbose, land, forced, method, directed, split, seed, samples = parseargs()
distances = zeros(length(vweights))
init_edges = Array{Int,2}(undef,0,0)
init_vweights = Vector{Float64}()
init_eweights = Vector{Float64}()
init_embed = Array{Float64,2}(undef,0,0)
v_to_l = Int[]
if land != -1
    init_edges = copy(edges)
    init_vweights = copy(vweights)
    init_eweights = copy(weights)
    init_embed = copy(embed)
    distances, embed, comm, edges, weights, vweights, v_to_l = landmarks(edges,
            weights, vweights, clusters, comm, embed, verbose, land, forced, method, directed)
end
if directed
    results = wGCL_directed(edges, weights, comm, embed, distances, vweights, init_vweights,
                v_to_l, init_edges, init_eweights, init_embed, split, seed, samples, verbose)
else
    results = wGCL(edges, weights, comm, embed, distances, vweights, init_vweights,
                v_to_l, init_edges, init_eweights, init_embed, split, seed, samples, verbose)
end
println(results)
