using Pkg
using ABCDGraphGenerator
using Random

for ξ in [0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65], ρ in [0.0], η = [1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0], d in [2], seed in 1:10
#for ξ in [0.25], η = [1.5], d in [2], ρ in [0.0], seed in 1:1 # disable this and enable above for larger sweep
    @show ξ, η, ρ, seed
    Random.seed!(seed)
    nout = 0
    n = 1000
    # in what follows n is number of non-outlier nodes
    n = n - nout

    τ₁ = 2.5
    d_min = 5 
    d_max = 50
    d_max_iter = 1000
    @info "Expected value of degree: $(ABCDGraphGenerator.get_ev(τ₁, d_min, d_max))"
    degs = ABCDGraphGenerator.sample_degrees(τ₁, d_min, d_max, n + nout, d_max_iter)
    open(io -> foreach(d -> println(io, d), degs), "degreefile_$(ξ)_$(η).txt", "w")

    τ₂ = 1.5
    c_min = 50
    c_max = 200
    c_max_iter = 1000
    @info "Expected value of community size: $(ABCDGraphGenerator.get_ev(τ₂, c_min, c_max))"
    coms = ABCDGraphGenerator.sample_communities(τ₂, ceil(Int, c_min / η), floor(Int, c_max / η), n, c_max_iter)
    @assert sum(coms) == n
    pushfirst!(coms, nout)

    p = ABCDGraphGenerator.ABCDParams(degs, coms, ξ, η, d, ρ)
    edges, clusters = ABCDGraphGenerator.gen_graph(p)
    open("networkfile_$(ξ)_$(η)_$(seed).txt", "w") do io
        for (a, b) in sort!(collect(edges))
            println(io, a, "\t", b)
        end
    end
    open("communityfile_$(ξ)_$(η)_$(seed).txt", "w") do io
        for (i, c) in enumerate(clusters)
            println(io, i, "\t", c)
        end
    end

    open("communitysizesfile_$(ξ)_$(η).txt", "w") do io
        comm_count = zeros(Int, length(coms))
        for c in clusters
            @assert length(c) > 0
            if 1 in c
                @assert length(c) == 1
            end
            for v in c
                comm_count[v] += 1
            end
        end
        println("eta is $η and empirically we have scaling of: ", extrema(comm_count[2:end] ./ coms[2:end]))
        foreach(d -> println(io, d), comm_count)
    end
end
