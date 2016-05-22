using Iterators
t = 10000
n = 1000
l = 50
lhat = 50
theta_a = 1/50
theta_v = 1/50
sa = linspace(0,3,7)
sv = linspace(0,3,7)
using Distributions
gamma_a = [1, 2]
gamma_b = [1/2, 1/3]
parameters = collect(product(t, n, l, lhat, theta_a, theta_v, sa, sv, gamma_a, gamma_b))

addprocs(24, topology=:master_slave)
addprocs([("yards.bio.upenn.edu", 25)], dir="/home/jakub/", exename="/home/jakub/bin/julia", topology=:master_slave)

@everywhere include("/home/jakub/Dropbox/coevolution/julia_sim/Coevolution.jl")

@everywhere function init_and_evolve(params)
    t, n, l, lhat, theta_a, theta_v, sa, sv, gamma_a, gamma_b = params
    dist = Gamma(gamma_a, gamma_b)
    kappa = rand(dist, l)
    A1 = Lineage(name = "antibody", n = n, s = sa, kappa = kappa,
    shat = sa, kappahat = rand(dist, lhat), randclone = randclone_n)
    A = Population(lineages = [A1], theta = theta_a)

    V1 = Lineage(name = "virus", n = n, s = -sv, kappa = kappa, kappahat = ones(0))
    V = Population(lineages = [V1], theta = theta_v)

    #evolve!(A, V, untilfixation = true )
    output = evolve!(A, V, tmax = t)

    output[:sa] = sa
    output[:sv] = sv
    output[:gamma_a] = gamma_a
    output[:gamma_b] = gamma_b
    output = output[output[:t] .> 2/theta_a,:]
    output = aggregate(output, [:sa, :sv, :name, :gamma_a, :gamma_b], mean)
    return output
end
@time result_stationary = reduce(vcat, pmap(init_and_evolve, parameters))
rmprocs(workers())
#end

# save data using R
using RCall
@rput(result_stationary)

R"save(result_stationary, file=\"figure_random_kappas.rdata\")"
# load data using R
#R"load(\"figure_stationary_random_kappas.rdata\")"
#@rget(result_stationary)
