using Iterators
t = 1000000
n = 1000
l = 50
lhat = 50
theta_a = 1/50
theta_v = 1/50
sa = linspace(0,3,7)
sv = linspace(0,3,7)

parameters = collect(product(t, n, l, lhat, theta_a, theta_v, sa, sv))

addprocs(24, topology=:master_slave)
addprocs([("yards.bio.upenn.edu", 25)], dir="/home/jakub/", exename="/home/jakub/bin/julia", topology=:master_slave)

@everywhere include("/home/jakub/Dropbox/coevolution/julia_sim/Coevolution.jl")

@everywhere function init_and_evolve(params)
    t, n, l, lhat, theta_a, theta_v, sa, sv = params

    A1 = Lineage(name = "antibody", n = n, s = sa, kappa = ones(l),
    shat = sa, kappahat = ones(lhat), randclone = randclone_n)
    A = Population(lineages = [A1], theta = theta_a)

    V1 = Lineage(name = "virus", n = n, s = -sv, shat = -sv,
      kappa = ones(l), kappahat = ones(lhat))
    V = Population(lineages = [V1], theta = theta_v)

    #evolve!(A, V, untilfixation = true )
    output = evolve!(A, V, tmax = t)

    output[:sa] = sa
    output[:sv] = sv

    output = output[output[:t] .> 2/theta_a,:]
    output = aggregate(output, [:sa, :sv, :name], mean)
    # output = @linq output |> where(:t .> 2/theta_a) |>
    #   by([:sa, :sv, :name], E = mean(:E), Ehat = mean(:Ehat),
    #    M2 = mean(:M2), M2hat = mean(:M2hat),
    #    F = mean(:F) , F2 = mean(:F2), EEcov = mean(:EEcov))
    return output
end
@time result_stationary = reduce(vcat, pmap(init_and_evolve, parameters))
rmprocs(workers())
#end

# save data using R
using RCall
@rput(result_stationary)
#params = DataFrame(t = t, n = n, l = l, lhat = lhat, theta_a = theta_a, theta_v = theta_v)
#@rput(params)
R"save(result_stationary, file=\"figure_moments.rdata\")"
