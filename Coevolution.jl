#module Coevolution
#export evolve!
#export Population

using Distributions
using DataFrames

type Clone
  g::Vector{Int8}             # genotype a vector of spins
  e::Float64                  # mean interaction of genotype with other population
  ghat::Vector{Int8}             # genotype a vector of spins
  ehat::Float64                  # mean interaction of genotype with other population
  n::Int64        # number of clones with this genotype
  f::Float64      # fitness of this genotype
  w::Float64  # weight for probability to be selected
  v::Float64  # interaction variance for e+ehat
  function Clone(g, ghat)
    new(g, 0.0, ghat, 0.0, 1, 0.0, 0.0, 0.0)
  end
end

function Base.copy(clone::Clone)
  newclone = Clone(copy(clone.g), copy(clone.ghat))
  #newclone.e = clone.e
  #newclone.ehat = clone.ehat
  newclone.n = clone.n
  #newclone.f = clone.f
  return newclone
end

function randclone_1(n::Int64, l, lhat) # makes a population if size n of 1 random genotypes of length l
  clones = Clone[]
  c = Clone(rand(-1:2:1, l), rand(-1:2:1, lhat))
  c.n = n
  push!(clones, c)
  return clones
end

function randclone_n(n::Int64, l, lhat) # makes a population of size n with n random genotypes of length l
  clones = Clone[]
  sizehint!(clones, n)
  for i in 1:n
    push!(clones, Clone(rand(-1:2:1, l), rand(-1:2:1, lhat)))
  end
  return clones
end

type Lineage
  clones::Vector{Clone}      # vector of clonal populations
  s::Float64                        # UNSCALED selection coefficient
  shat::Float64 # UNSCALED selection coefficient for conserved region
  kappa::Vector{Float64}  # site specific effects
  kappahat::Vector{Float64}  # site specific effects
  e02::Float64
  e02hat::Float64
  name::ASCIIString
  beta0::Float64 # nonlinear fitness parameters
  logM::Float64
  Estar0::Float64
  function Lineage(;n = 500, s = 0.0, kappa = Float64[],
                    shat = 0.0, kappahat = Float64[],
                    randclone = randclone_n,
                    name = "lineage", beta = 0.0, logM = 0, Estar = 0)
    clones = randclone(n, length(kappa), length(kappahat))
    e02 = sum(kappa.^2)
    if e02 > 0.0
      s_unscaled = s/sqrt(e02) # scaled by n later
    else
      s_unscaled = 0.0
      e02 = 1 # this helps with printing normalized values
    end
    e02hat = sum(kappahat.^2)
    if e02hat > 0.0
      shat_unscaled = shat/sqrt(e02hat)
    else
      shat_unscaled = 0.0
      e02hat = 1
    end
    beta0 = beta/sqrt(e02+e02hat)
    Estar0 = Estar*sqrt(e02+e02hat)
    new(clones, s_unscaled, shat_unscaled, kappa, kappahat, e02, e02hat, name, beta0, logM, Estar0)
  end
end


function Base.copy(lin::Lineage)
  newlin = Lineage(n=0)
  newlin.clones = Clone[]
  for clone in lin.clones
    push!(newlin.clones, copy(clone))
  end
  newlin.s = lin.s
  newlin.shat = lin.shat
  newlin.kappa = lin.kappa
  newlin.kappahat = lin.kappahat
  newlin.e02 = lin.e02
  newlin.e02hat = lin.e02hat
  newlin.name = lin.name
  return newlin
end


type Population
  lineages::Vector{Lineage}      # vector of clonal populations
  n::Int64                          # total population size
  f::Float64                        # mean fitness
  l::Int64                          # number of sites
  lhat::Int64  # length of conserved region
  theta_pois::Distributions.Poisson # poisson random number generator for mutations
  theta_pois_hat::Distributions.Poisson # poisson random number generator for mutations
  function Population(;theta = 1/30, lineages = [])
    ls = unique([(length(lin.kappa),length(lin.kappahat)) for lin in lineages])
    if length(ls) > 1
      error("kappas in all lineages should have same length")
    end
    l, lhat = ls[1]
    theta_pois = Poisson(theta*l)
    theta_pois_hat = Poisson(theta*lhat)
    # count total pop size from lineages
    n = 0
    for lin in lineages
      for clone in lin.clones
        n += clone.n
      end
    end
    # rescale s in each lineage
    # for lin in lineages
    #   lin.s = lin.s/n
    #   lin.shat = lin.shat/n
    # end

    new(lineages, n, 0.0,  l, lhat, theta_pois, theta_pois_hat)
  end
end

function Base.copy(pop::Population)
  newpop = Population(lineages = pop.lineages)
  newpop.theta_pois = pop.theta_pois
  newpop.theta_pois_hat = pop.theta_pois_hat
  newpop.lineages = Lineage[]
  for lin in pop.lineages
    push!(newpop.lineages, copy(lin))
  end
  return newpop
end

function interaction(gA::Vector{Int8}, gV::Vector{Int8}, kappa::Vector{Float64}, l::Int64)
  e = 0.0
  @simd for i = 1:l
    @fastmath @inbounds e += kappa[i]*gA[i]*gV[i]
    #e += kappa[i]*gA[i]*gV[i]
  end
  return e
end

function interaction(gA::Vector{Int8}, kappa::Vector{Float64}, l::Int64)
  e = 0.0
  @simd for i = 1:l
    @fastmath @inbounds e += kappa[i]*gA[i]
    #e += kappa[i]*gA[i]
  end
  return e
end

function interaction!(A::Population, V::Population) # calculates mean interactions between populations
  ehatmean = 0.0
  for lineage in A.lineages
    for clone in lineage.clones
      clone.e = 0.0
      #clone.ehat = sum(clone.ghat .* A.kappahat)
      clone.ehat = interaction(clone.ghat, lineage.kappahat, A.lhat)
      ehatmean += clone.ehat*clone.n/A.n
    end
  end
  for lineage in V.lineages
    for clone in lineage.clones
      clone.e = 0.0
      clone.ehat = ehatmean
    end
  end
  for Alineage in A.lineages, clonea in Alineage.clones
    for Vlineage in V.lineages, clonev in Vlineage.clones
      #e = sum(clonea.g .* clonev.g .* A.kappa)
      e = interaction(clonea.g, clonev.g, Alineage.kappa, A.l)
      clonea.e += e*clonev.n/V.n
      clonev.e += e*clonea.n/A.n
    end
  end
  # for lineage in A.lineages
  #   for clone in lineage.clones
  #     clone.e /= V.n
  #   end
  # end
  # for lineage in V.lineages
  #   for clone in lineage.clones
  #     clone.e /= A.n
  #   end
  # end

end

function fitness_linear!(pop::Population)
  pop.f = 0.0
  for lineage in pop.lineages
    for clone in lineage.clones
      clone.f = lineage.s/pop.n*clone.e + lineage.shat/pop.n*clone.ehat
      pop.f += clone.f*clone.n/pop.n
    end
  end
end

function fitness_linear!(A::Population, V::Population)
  interaction!(A, V)
  fitness_linear!(A), fitness_linear!(V)
end

function interaction_var!(A::Population, V::Population) # calculates mean interactions between populations
  for lineage in A.lineages
    for clone in lineage.clones
      clone.v = 0.0
    end
  end
  for lineage in V.lineages
    for clone in lineage.clones
      clone.v = 0.0
    end
  end
  for Alineage in A.lineages, clonea in Alineage.clones
    ehat = interaction(clonea.ghat, Alineage.kappahat, A.lhat)
    for Vlineage in V.lineages, clonev in Vlineage.clones
      e = interaction(clonea.g, clonev.g, Alineage.kappa, A.l)
      clonea.v += (e - clonea.e)^2*clonev.n/V.n
      clonev.v += (e + ehat - clonev.e - clonev.ehat)^2*clonea.n/A.n
    end
  end
end

function fitness_nonlinear!(A::Population, V::Population)
  interaction!(A, V)
  interaction_var!(A, V)
  A.f = 0.0
  for lineage in A.lineages
    for clone in lineage.clones
      clone.f = -log(1+exp(-lineage.beta0*(clone.e + clone.ehat + sqrt(2*lineage.logM*clone.v) - lineage.Estar0)))*lineage.s/A.n/lineage.beta0*(1+exp(-lineage.beta0*lineage.Estar0))
      A.f += clone.f*clone.n/A.n
    end
  end
  V.f = 0.0
  for lineage in V.lineages
    for clone in lineage.clones
      clone.f = -log(1+exp(-lineage.beta0*(clone.e + clone.ehat - lineage.Estar0)))*lineage.s/V.n/lineage.beta0*(1+exp(-lineage.beta0*lineage.Estar0))
      V.f += clone.f*clone.n/V.n
    end
  end
end

function selection!(pop::Population)
  Z = 0.0
  i = 0
  for lin in pop.lineages
    for clone in lin.clones
      clone.w = clone.n*exp(clone.f - pop.f)
      Z += clone.w
      i += 1
      #println("$(lin.name) $(clone.n)")
    end
  end
  nr = pop.n
  deadlini = Int64[]
  for (lini, lin) in enumerate(pop.lineages)
    dead = Int64[]
    for (ci, clone) in enumerate(lin.clones)
      if i == 1
        clone.n = nr
      else
        clone.n = rand(Binomial(nr, clone.w/Z))
        Z -= clone.w
        nr -= clone.n
        i -= 1
      end
      if clone.n == 0
        push!(dead, ci)
      end
    end
    deleteat!(lin.clones, dead)
    if length(lin.clones) == 0
      push!(deadlini, lini)
    end
  end
  deleteat!(pop.lineages, deadlini)
  return length(pop.lineages) == 1
end

# function selection!(pop::Population)
#   p = Float64[]
#   for lin in pop.lineages
#     for clone in lin.clones
#       push!(p, clone.n*exp(clone.f - pop.f))
#     end
#   end
#   n = rand(Multinomial(pop.n, p/sum(p)))
#   ni = 1
#   deadlini = Int64[]
#   for (lini, lin) in enumerate(pop.lineages)
#     dead = Int64[]
#     for (ci, clone) in enumerate(lin.clones)
#       clone.n = n[ni]
#       ni += 1
#       if clone.n == 0
#           push!(dead, ci)
#       end
#     end
#     deleteat!(lin.clones, dead)
#     if length(lin.clones) == 0
#       push!(deadlini, lini)
#     end
#   end
#   # deadlin = pop.lineages[deadlini]
#   deleteat!(pop.lineages, deadlini)
#   return length(pop.lineages) == 1
# end


function random_clone(lineages::Vector{Lineage}, n::Int64)
  k = rand(1:n)
  x = 0
  for lin in lineages
    for (i, clone) in enumerate(lin.clones)
      x += clone.n
      if k <= x
        return (lin, i, clone)
      end
    end
  end
end

function mutate!(pop::Population)
  if pop.l > 0
    for k = 1:rand(pop.theta_pois)
      lineage, clone_i, clone = random_clone(pop.lineages, pop.n)
      clone.n -= 1
      mutated_g = copy(clone.g)
      i = rand(1:pop.l)
      mutated_g[i] = -mutated_g[i]
      newclone = Clone(mutated_g, clone.ghat)
      if clone.n == 0
        deleteat!(lineage.clones, clone_i)
      end
      push!(lineage.clones, newclone)
    end
  end
  if pop.lhat > 0
    for k = 1:rand(pop.theta_pois_hat)
      lineage, clone_i, clone = random_clone(pop.lineages, pop.n)
      clone.n -= 1
      mutated_g = copy(clone.ghat)
      i = rand(1:pop.lhat)
      mutated_g[i] = -mutated_g[i]
      newclone = Clone(clone.g, mutated_g)
      if clone.n == 0
        deleteat!(lineage.clones, clone_i)
      end
      push!(lineage.clones, newclone)
    end
  end
end

function linmoment(lin::Lineage, popn::Int64, xs::Symbol, m::Float64 = 0.0, p::Int64 = 1)
  x = 0.0
  for clone in lin.clones
    x += clone.n/popn*(getfield(clone, xs) - m)^p
  end
  return x
end

function ecov(lin::Lineage, popn::Int64, e::Float64 = 0.0, ehat::Float64 = 0.0)
  x = 0.0
  for clone in lin.clones
    x += clone.n/popn*(clone.e - e)*(clone.ehat - ehat)
  end
  return x
end

function summary!(output, pop::Population, t, scaled = true)
  Fpop = 0.0
  Epop = 0.0
  Ehatpop = 0.0
  for lin in pop.lineages
    Fpop += linmoment(lin, pop.n, :f)
    Epop += linmoment(lin, pop.n, :e)
    Ehatpop += linmoment(lin, pop.n, :ehat)
  end
  for lin in pop.lineages
    x     = linmoment(lin, pop.n, :n, 0.0, 0)
    F     = linmoment(lin, pop.n, :f)/x
    F2    = linmoment(lin, pop.n, :f, Fpop, 2)/x
    E     = linmoment(lin, pop.n, :e)/x
    Ehat  = linmoment(lin, pop.n, :ehat)/x
    M2    = linmoment(lin, pop.n, :e, Epop, 2)/x
    M2hat = linmoment(lin, pop.n, :ehat, Ehatpop, 2)/x
    M3    = linmoment(lin, pop.n, :e, Epop, 3)/x
    M3hat = linmoment(lin, pop.n, :ehat, Ehatpop, 3)/x
    M4    = linmoment(lin, pop.n, :e, Epop, 4)/x
    M4hat = linmoment(lin, pop.n, :ehat, Ehatpop, 4)/x
    EEcov = ecov(lin, pop.n, Epop, Ehatpop)/x
    if scaled then
      d = @data([t, lin.name, x, F*pop.n, F2*pop.n^2, E/sqrt(lin.e02), M2/lin.e02,
        Ehat/sqrt(lin.e02hat), M2hat/lin.e02hat,
        EEcov/sqrt(lin.e02)/sqrt(lin.e02hat),
        M3*lin.e02^(-3/2), M3hat*lin.e02hat^(-3/2), M4*lin.e02^(-2), M4hat*lin.e02hat^(-2)])
    else
      d = @data([t, lin.name, x, F*pop.n, F2*pop.n^2, E, M2,
        Ehat, M2hat, EEcov, M3, M3hat, M4, M4hat])
    end
    push!(output, d)
  end
end

function evolve!(A, V; tmax = 50, untilfixation = false, fitness! = fitness_linear!, scaled = true, tout = A.n)
  # check lengths are consistent betwewen A and V
  if (A.l != V.l) || (A.lhat != A.lhat)
    error("genotype lengths of populations are not equal")
  end
  output = DataFrame(t = Float64[], name = ASCIIString[], frequency = Float64[],
                     F = Float64[], F2 = Float64[],
                     E = Float64[], M2 = Float64[],
                     Ehat = Float64[], M2hat = Float64[], EEcov = Float64[],
                     M3 = Float64[], M3hat = Float64[], M4 = Float64[], M4hat = Float64[])
  fitness!(A, V)
  t = 0
  summary!(output, A, t/A.n, scaled), summary!(output, V, t/A.n, scaled)
  fixed = false
  tmax = tmax*A.n
  while untilfixation || t < tmax
    for tn in 1:tout
      t += 1
      fixed = selection!(A)
      selection!(V)
      mutate!(A), mutate!(V)
      fitness!(A, V)
      if untilfixation && fixed
        break
      end
    end
    if !untilfixation
      summary!(output, A, t/A.n, scaled), summary!(output, V, t/A.n, scaled)
      #println(t/A.n)
    elseif fixed
      summary!(output, A, t/A.n, scaled), summary!(output, V, t/A.n, scaled)
      #println(t/A.n)
      break
    end
      # push!(output, @data([t,
      #   A.n*A.f, V.n*V.f,
      #   E/sqrt(A.e02), Ma2/A.e02, Mv2/V.e02,
      #   Ehat/sqrt(A.e02hat), Ma2hat/A.e02hat, Mv2hat/A.e02hat]))
  end
  return output
end


#end
