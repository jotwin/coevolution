n = 1000
l = 50
lhat = 50
theta_a = 1/50
theta_v = 1/50

addprocs(24, topology=:master_slave)
addprocs([("yards.bio.upenn.edu", 25)], dir="/home/jakub/", exename="/home/jakub/bin/julia", topology=:master_slave)

@everywhere include("/home/jakub/Dropbox/coevolution/julia_sim/Coevolution.jl")

@everywhere function init_and_evolve(parameters)
  n, l, lhat, theta_a, theta_v, sa, sv = parameters
  A1 = Lineage(name = "antibody", n = n, s = sa, kappa = ones(l),
    shat = sa, kappahat = ones(lhat), randclone = randclone_1)
  A = Population(lineages = [A1], theta = theta_a)

  V1 = Lineage(name = "virus", n = n, s = -sv, kappa = ones(l), kappahat = ones(0), randclone = randclone_1)
  V = Population(lineages = [V1], theta = theta_v)

  evolve!(A, V, tmax = 1)
  output = DataFrame(t = Float64[], name = ASCIIString[], frequency = Float64[],
                     F = Float64[], F2 = Float64[],
                     E = Float64[], M2 = Float64[],
                     Ehat = Float64[], M2hat = Float64[], EEcov = Float64[],
                     M3 = Float64[], M3hat = Float64[], M4 = Float64[], M4hat = Float64[])
  fitness! = fitness_linear!
  fitness!(A, V)
  summary!(output, A, 0)
  summary!(output, V, 0)
  A0 = copy(A)
  V0 = copy(V)
  Ahistory = Population[]
  Vhistory = Population[]
  tmax = convert(Int64, 2/theta_a)
  for t = 1:tmax
    for tn in 1:A.n
      selection!(A), selection!(V)
      mutate!(A), mutate!(V)
      fitness!(A, V)
    end
    push!(Ahistory, copy(A))
    push!(Vhistory, copy(V))
  end
  for t = 1:tmax
    Vh = Vhistory[t]
    fitness!(A0, Vh)
    summary!(output, A0, -t)
    summary!(output, Vh, -t)
  end
  for t = 1:tmax
    Ah = Ahistory[t]
    fitness!(Ah, V0)
    summary!(output, V0, t)
    summary!(output, Ah, t)
  end
  #output = aggregate(output, [:name, :t], mean)
  output[:theta_a] = theta_a
  output[:theta_v] = theta_v
  output[:sa] = sa
  output[:sv] = sv
  return output
end

parameters1 = rep((n, l, lhat, theta_a, theta_v, 1, 2), 1000)
@time output1 = reduce(vcat, pmap(init_and_evolve, parameters1))
parameters2 = rep((n, l, lhat, theta_a, theta_v, 2, 1), 1000)
@time output2 = reduce(vcat, pmap(init_and_evolve, parameters2))

#init_and_evolve(parameters1[1])

# output1 = init_and_evolve(n, l, lhat, theta_a, theta_v, 1, 2)
# output2 = init_and_evolve(n, l, lhat, theta_a, theta_v, 2, 1)
output = vcat(output1, output2)

output = aggregate(output, [:name, :sa, :sv, :theta_a, :theta_v, :t], mean)
# save data using R
using RCall
@rput(output)
R"save(output, file=\"figure_timeshift_nonstationary.rdata\")"
# load data using R
# R"load(\"figure_timeshift_nonstationary.rdata\")"
# @rget(output)
#
# using DataFrames
# using DataFramesMeta
# using Gadfly
# import Lazy
# theme  = Theme(major_label_font_size = 9pt,
#   minor_label_font_size = 6pt,
#   key_label_font_size = 6pt,
#   key_title_font_size = 9pt,
#   default_point_size = 1pt,
#   highlight_width = 0pt)



### data wrangling
# outputg = Lazy.@> begin
#   output
#   @where(:name .== "virus", t .> 0.5)
#   @transform(t = -:t)
#   vcat(@where(output, :name .== "antibody"))
#   @transform(Fa = :sa*:E_mean)
# end
#
# p = plot(outputg,
#   layer(x = :t, y= :Fa, Geom.line), theme)
# draw(PDF("figures/timeshift_nonstationary.pdf", 8.3cm, 6cm), p)
