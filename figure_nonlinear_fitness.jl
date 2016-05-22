using Iterators
t = 100000
n = 1000
l = 50
lhat = 0
theta_a = 1/50
theta_v = 1/50
sa = linspace(0,3,7)
sv = 1
beta = [.1, 1, 10]
logM = log([1, 10, 100, 1000])
Estar = [0, 5]

parameters = collect(product(t, n, l, lhat, theta_a, theta_v, sa, sv, beta, logM, Estar))

addprocs(24, topology=:master_slave)
addprocs([("yards.bio.upenn.edu", 25)], dir="/home/jakub/", exename="/home/jakub/bin/julia", topology=:master_slave)

@everywhere include("/home/jakub/Dropbox/coevolution/julia_sim/Coevolution.jl")

@everywhere function init_and_evolve(params)
    t, n, l, lhat, theta_a, theta_v, sa, sv, beta, logM, Estar = params

    A1 = Lineage(name = "antibody", n = n, s = sa, kappa = ones(l),
    shat = sa, kappahat = ones(lhat), randclone = randclone_n,
    beta = beta, logM = logM, Estar = Estar)
    A = Population(lineages = [A1], theta = theta_a)

    V1 = Lineage(name = "virus", n = n, s = -sv, kappa = ones(l), kappahat = ones(0),
    beta = beta, logM = logM, Estar = Estar)
    V = Population(lineages = [V1], theta = theta_v)

    output = evolve!(A, V, tmax = t, fitness! = fitness_nonlinear!)

    output[:sa] = sa
    output[:sv] = sv
    output[:beta] = beta
    output[:logM] = logM
    output[:Estar] = Estar
    output = output[output[:t] .> 2/theta_a,:]
    output = aggregate(output, [:sa, :sv, :name], mean)
    return output
end
@time result_stationary = reduce(vcat, pmap(init_and_evolve, parameters))
rmprocs(workers())
#end

# save data using R
using RCall
@rput(result_stationary)
R"save(result_stationary, file=\"figure_nonlinear_fitness.rdata\")"
