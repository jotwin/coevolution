t = 10000
n = 1000
l = 50
lhat = 50
theta_a = 1/50
theta_v = 1/50
sa = 1
sv = 1

include("Coevolution.jl")

A1 = Lineage(name = "antibody", n = n, s = sa, kappa = ones(l),
  shat = sa, kappahat = ones(lhat), randclone = randclone_n)
A = Population(lineages = [A1], theta = theta_a)

V1 = Lineage(name = "virus", n = n, s = -sv, kappa = ones(l),
  shat = -sv, kappahat = ones(lhat))
V = Population(lineages = [V1], theta = theta_v)

output = evolve!(A, V, tmax = t)

# save data using R
using RCall
@rput(output)
R"save(output, file=\"figure_stationary_timecorrelation.rdata\")"
