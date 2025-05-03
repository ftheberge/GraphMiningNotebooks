using Pkg
using ABCDGraphGenerator
using Random

@info "Usage: julia abcd_sampler.jl config_filename"
@info "For the syntax of config_filename see example_config.toml file"

filename = ARGS[1]
conf = Pkg.TOML.parsefile(filename)
isempty(conf["seed"]) || Random.seed!(parse(Int, conf["seed"]))

nout = haskey(conf, "nout") ? parse(Int, conf["nout"]) : 0

μ = haskey(conf, "mu") ? parse(Float64, conf["mu"]) : nothing
ξ = haskey(conf, "xi") ? parse(Float64, conf["xi"]) : nothing
if !(isnothing(μ) || isnothing(ξ))
    throw(ArgumentError("inconsistent data: only μ or ξ may be provided"))
end

if !isnothing(μ) && nout > 0
    throw(ArgumentError("μ is not supported with outliers"))
end

n = parse(Int, conf["n"])

if nout > n
    throw(ArgumentError("number of outliers cannot be larger than graph size"))
end

islocal = haskey(conf, "islocal") ? parse(Bool, conf["islocal"]) : false
if islocal && nout > 0
    throw(ArgumentError("local graph is not supported with outliers"))
end

isCL = parse(Bool, conf["isCL"])
if isCL && nout > 0
    throw(ArgumentError("Chung-Lu graph is not supported with outliers"))
end

# in what follows n is number of non-outlier nodes
n = n - nout

τ₁ = parse(Float64, conf["t1"])
d_min = parse(Int, conf["d_min"])
d_max = parse(Int, conf["d_max"])
d_max_iter = parse(Int, conf["d_max_iter"])
@info "Expected value of degree: $(ABCDGraphGenerator.get_ev(τ₁, d_min, d_max))"
degs = ABCDGraphGenerator.sample_degrees(τ₁, d_min, d_max, n + nout, d_max_iter)
open(io -> foreach(d -> println(io, d), degs), conf["degreefile"], "w")

τ₂ = parse(Float64, conf["t2"])
c_min = parse(Int, conf["c_min"])
c_max = parse(Int, conf["c_max"])
c_max_iter = parse(Int, conf["c_max_iter"])
@info "Expected value of community size: $(ABCDGraphGenerator.get_ev(τ₂, c_min, c_max))"
coms = ABCDGraphGenerator.sample_communities(τ₂, c_min, c_max, n, c_max_iter)
if nout > 0
    pushfirst!(coms, nout)
end
open(io -> foreach(d -> println(io, d), coms), conf["communitysizesfile"], "w")

p = ABCDGraphGenerator.ABCDParams(degs, coms, μ, ξ, isCL, islocal, nout > 0)
edges, clusters = ABCDGraphGenerator.gen_graph(p)
open(conf["networkfile"], "w") do io
    for (a, b) in sort!(collect(edges))
        println(io, a, "\t", b)
    end
end
open(conf["communityfile"], "w") do io
    for (i, c) in enumerate(clusters)
        println(io, i, "\t", c)
    end
end
