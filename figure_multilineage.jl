
using Distributions

begin
  lineages = []
  for i = 1:10
    lin = A1 = Lineage(name = "antibody $(i)", n = 1000, s = 1, kappa = rand(Gamma(1,1/3), 50),
          shat = 1, kappahat = rand(Gamma(1,1/3), 50), randclone = randclone_n)
    push!(lineages, lin)
  end

  A = Population(lineages = lineages, theta = 1/50)

  V1 = Lineage(name = "virus", n = 10000, s = -1, kappa = ones(50), kappahat = ones(0))
  V = Population(lineages = [V1], theta = 1/50)

  #evolve!(A, V, untilfixation = true )
  @time output = evolve!(A, V, tmax = .5, tout = 10)
end
using RCall
@rput(output)
R"library(ggplot2)"
R"library(dplyr)"
R"qplot(data=output %>% filter(name!=\"virus\"), x=t, y=frequency, geom=\"area\" , fill=name)"
R"save(output, file=\"figure_multilineage.rdata\")"
R"ggsave(\"figures/multilineage.pdf\", height=4, width=5)"
