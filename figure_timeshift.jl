#do
using Iterators
t = 10000
n = 1000
l = 50
lhat = 50
theta_a = [1/50, 3/50]
theta_v = [1/50, 3/50]
sa = [0, 2]
sv = [0, 2]

parameters = collect(product(t, n, l, lhat, theta_a, theta_v, sa, sv))

addprocs(24, topology=:master_slave)
addprocs([("yards.bio.upenn.edu", 24)], dir="/home/jakub/", exename="/home/jakub/bin/julia", topology=:master_slave)
# this produces a lot of warnings
@everywhere include("/home/jakub/Dropbox/coevolution/julia_sim/Coevolution.jl")


@everywhere function evolve_time_shift!(A, V, tmax, dtmax, tcut, fitness! = fitness_linear!)
  # check lengths are consistent betwewen A and V
  if A.l != V.l
    error("genotype lengths of populations are not equal")
  end
  output = DataFrame(t = Float64[], name = ASCIIString[], frequency = Float64[],
                     F = Float64[], F2 = Float64[],
                     E = Float64[], M2 = Float64[],
                     Ehat = Float64[], M2hat = Float64[], EEcov = Float64[],
                     M3 = Float64[], M3hat = Float64[], M4 = Float64[], M4hat = Float64[])
  Ahistory = Population[]
  Vhistory = Population[]
  fitness!(A, V)
  for t = 1:tmax
    for tn in 1:A.n
      selection!(A), selection!(V)
      mutate!(A), mutate!(V)
      fitness!(A, V)
    end
    #println(t)
    #summary!(output, A, t), summary!(output, V, t)
    push!(Ahistory, copy(A))
    push!(Vhistory, copy(V))
  end
  for tA = (tcut+dtmax):(tmax-dtmax)
    for tV = (tA-dtmax):(tA+dtmax)
      Ah = Ahistory[tA]
      Vh = Vhistory[tV]
      fitness!(Ah, Vh)
      summary!(output, Ah, tV-tA), summary!(output, Vh, tA-tV)
    end
  end
  return output
end

@everywhere function init_and_evolve(params)
    t, n, l, lhat, theta_a, theta_v, sa, sv = params

    A1 = Lineage(name = "antibody", n = n, s = sa, kappa = ones(l),
    shat = sa, kappahat = ones(lhat), randclone = randclone_n)
    A = Population(lineages = [A1], theta = theta_a)

    V1 = Lineage(name = "virus", n = n, s = -sv, kappa = ones(l), kappahat = ones(0))
    V = Population(lineages = [V1], theta = theta_v)

    #evolve!(A, V, untilfixation = true )
    output = evolve_time_shift!(A, V, t, 100, 100)

    output[:sa] = sa
    output[:sv] = sv
    output[:lhat] = lhat
    output[:l] = l
    output[:theta_a] = theta_a
    output[:theta_v] = theta_v
    output = aggregate(output, [:sa, :sv, :name, :lhat, :l, :theta_a, :theta_v, :t], mean)
    return output
end
@time result_stationary = reduce(vcat, pmap(init_and_evolve, parameters))
rmprocs(workers())

#init_and_evolve(parameters[])

# save data using R
using RCall
@rput(result_stationary)
R"save(result_stationary, file=\"figure_timeshift.rdata\")"
# load data using R
# R"load(\"figure_timeshift.rdata\")"
# @rget(result_stationary)
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
#
#
#
#
# ### data wrangling, put M2A M2V as new column
# result_stationary_m2 = Lazy.@> begin
#   result_stationary
#   @select(:name, :sa, :sv, :theta_a, :theta_v, :lhat, :t, :M2_mean)
#   unstack(:name, :M2_mean)
#   @select(:sa, :sv, :theta_a, :theta_v, :lhat, :t, M2A = :antibody, M2V = :virus)
#   join(@where(result_stationary, :name .== "antibody"), on = [:sa, :sv, :theta_a, :theta_v, :lhat, :t])
# end
# # # add conversed region
# # result_stationary_m2 = Lazy.@> begin
# #   result_stationary
# #   @select(:name, :sa, :sv, :theta_a, :theta_v, :lhat, :t, :M2hat_mean)
# #   unstack(:name, :M2hat_mean)
# #   @select(:sa, :sv, :theta_a, :theta_v, :lhat, :t, M2Ahat = :antibody, M2Vhat = :virus)
# #   join(result_stationary_m2, on = [:sa, :sv, :theta_a, :theta_v, :lhat, :t])
# # end
# # # stack regions
# # result_stationary_m2v = @linq result_stationary_m2 |>
# #   select(:sa, :sv, :theta_a, :theta_v, :lhat, :t, E = :E_mean, :M2A, :M2V) |>
# #   transform(region = "variable")
# # result_stationary_m2c = @linq result_stationary_m2 |>
# #   select(:sa, :sv, :theta_a, :theta_v, :lhat, :t, E = :Ehat_mean, M2A = :M2Ahat, M2V = :M2Vhat) |>
# #   transform(region = "conserved")
# # result_stationary_m2vc = vcat(result_stationary_m2v, result_stationary_m2c)
# # write theory for t > 0, t <0, theta_a != theta_v and theta_a == theta_v
# result_stationary1 = @linq result_stationary_m2 |>
#   where(:t .>= 0, :theta_a .!= :theta_v) |>
#   transform(E_theory = 1/2*:M2V.*:sv./(:theta_a-:theta_v).*exp(-2.*:theta_a.*:t) +
#     1/2*(:M2A.*:sa.*(:theta_a-:theta_v)-2*(:theta_a.*:M2V.*:sv))./(:theta_a.^2-:theta_v.^2).*exp(-2.*:theta_v.*:t))
# result_stationary2 = @linq result_stationary_m2 |>
#   where(:t .< 0, :theta_a .!= :theta_v) |>
#   transform(E_theory = 1/2*:M2A.*:sa./(:theta_a-:theta_v).*exp(2.*:theta_v.*:t) +
#     1/2*(:M2V.*:sv.*(:theta_v-:theta_a)-2*(:theta_v.*:M2A.*:sa))./(:theta_a.^2-:theta_v.^2).*exp(2.*:theta_a.*:t))
# result_stationary3 = @linq result_stationary_m2 |>
#   where(:t .>= 0, :theta_a .== :theta_v) |>
#   transform(E_theory = ((:M2A.*:sa-:M2V.*:sv)/4./:theta_a-:M2V.*:sv.*:t).*exp(-2*:theta_a.*:t))
# result_stationary4 = @linq result_stationary_m2 |>
#   where(:t .< 0, :theta_a .== :theta_v) |>
#   transform(E_theory = ((:M2A.*:sa-:M2V.*:sv)/4./:theta_a-:M2A.*:sa.*:t).*exp(2*:theta_a.*:t))
# result_stationary_theory = vcat(result_stationary1, result_stationary2, result_stationary3, result_stationary4)
# # plot E vs time-shift and theory
# p = plot(@where(result_stationary_theory, :sa .== 1, :sv .== 2, :theta_v .== 2/50, :lhat .== 50),
#   layer(x = :t, y= :E_theory, color = :theta_a, Geom.line),
#   layer(x = :t, y= :E_mean, color = :theta_a, Geom.point), theme)
# draw(PDF("figures/timeshift1.pdf", 8.3cm, 6cm), p)
# p = plot(@where(result_stationary_theory, :theta_v .== 2/50, :lhat .== 50),
#   Geom.subplot_grid(layer(x = :t, y= :E_mean, color = :theta_a, xgroup = :sa, ygroup = :sv, Geom.point),
#   layer(x = :t, y= :E_theory, color = :theta_a, xgroup = :sa, ygroup = :sv, Geom.line)),
#   theme)
# draw(PDF("figures/timeshift1.pdf", 8.0inch, 8.0inch), p)
#
#
#
# #### plot mean E at dt = 0 vs theory
# result_stationary2 = @linq result_stationary_m2 |>
#   where(:name .== "antibody", :t .== 0, :sa .== 2, :sv .== 1, :lhat .== 0) |>
#   transform(E = :E_mean, E_theory = (:sa.*:M2A-:sv.*:M2V)/2./(:theta_a+:theta_v))
#
# plot(result_stationary2, layer(x = "theta_v", y= "E", color = "theta_a", Geom.point),
#   layer(x = "theta_v", y= "E_theory", color = "theta_a", Geom.line))
