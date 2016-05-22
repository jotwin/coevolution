using Iterators

n = 1000
l = 50
lhat = 50
theta_a = [1/50, 2/50]
theta_v = [0.1/50, 1/50, 2/50, 3/50, 4/50, 5/50]
sa = [.5, .75, 1, 2]
sv = 1
id = 1:1000
bnab = [true, false]

p = collect(product(n, l, lhat, theta_a, theta_v, sa, sv, id, bnab))

addprocs(22, topology=:master_slave)
addprocs([("yards.bio.upenn.edu", 23)], dir="/home/jakub/", exename="/home/jakub/bin/julia", topology=:master_slave)

@everywhere include("/home/jakub/Dropbox/coevolution/julia_sim/Coevolution.jl")

@everywhere function LA1V1(alin::Lineage, A::Population, V::Population)
  fitness_linear!(A, V)
  emean = 0.0
  for lin in A.lineages
    emean += linmoment(lin, A.n, :e)
  end
  x = 0.0
  for vlin in V.lineages, vclone in vlin.clones
    for aclone in alin.clones
      e = interaction(aclone.g, vclone.g, alin.kappa, A.l)
      x += aclone.n/A.n*(e - emean)*vclone.n/V.n*(vclone.e - emean)
    end
  end
  return x
end

@everywhere function scale_clones!(clones, f)
  dead = []
  for (ci, clone) in enumerate(clones)
    clone.n = round(Int64, clone.n*f)
    if clone.n == 0
      push!(dead, ci)
    end
  end
  deleteat!(clones, dead)
end

@everywhere function init_and_evolve(parameters)
  n, l, lhat, theta_a, theta_v, sa, sv, id, bnab = parameters

  # initial evolution of resident and virus
  Alinres = Lineage(name = "antibody resident", n = n, s = sa, kappa = ones(l),
    shat = 0, kappahat = zeros(lhat), randclone = randclone_1)
  Ares = Population(lineages = [Alinres], theta = theta_a)
  V1 = Lineage(name = "virus", n = n, s = -sv, kappa = ones(l),
    shat = 0, kappahat = zeros(lhat), randclone = randclone_1)
  V = Population(lineages = [V1], theta = theta_v)
  Alininv = Lineage(name = "antibody invader", n = n, s = sa*(1-1*bnab), kappa = (1-1*bnab)*ones(l),
    shat = sa*bnab, kappahat = bnab*ones(lhat), randclone = randclone_1)
  Ainv = Population(lineages = [Alininv], theta = theta_a)

  for t = 1:floor(Int64, 50), tn in 1:n
    fitness_linear!(Ares, V)
    selection!(Ares)
    selection!(V)
    fitness_linear!(Ainv, V)
    selection!(Ainv)
    mutate!(Ares)
    mutate!(Ainv)
    mutate!(V)
  end

  # combine resident and invader
  freq0 = 0.10
  scale_clones!(Alinres.clones, 1-freq0)
  scale_clones!(Alininv.clones, freq0)

  A = Population(lineages = [Alinres, Alininv], theta = theta_a)
  la1v1 = LA1V1(Alininv, A, V)
  # evolve until fixation, record initial state
  A0 = copy(A)
  V0 = copy(V)
  output = evolve!(A, V, untilfixation = true, scaled = false)
  for r = 1:99
    A = copy(A0)
    V = copy(V0)
    output2 = evolve!(A, V, untilfixation = true, scaled = false)
    output = vcat(output, output2[output2[:t] .> 0, :])
  end
  output[:id] = id
  output[:LA1V1] = la1v1
  output[:sa] = sa
  output[:sv] = sv
  output[:theta_a] = theta_a
  output[:theta_v] = theta_v
  output[:bnab] = bnab
  return output
end

@time output = reduce(vcat, pmap(init_and_evolve, p))
rmprocs(workers())

#init_and_evolve(p[1])

using RCall
@rput(output)
R"save(output, file=\"/data/jakub/figure_fixation.rdata\")"
