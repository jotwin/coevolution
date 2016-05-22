# compete two lineages
A1 = Lineage(name = "antibody lin 1", n = 500, s = 1, kappa = ones(25),
            shat = 1, kappahat = ones(25), randclone = randclone_1)
A2 = Lineage(name = "antibody lin 2", n = 500, s = 1, kappa = ones(25),
                         shat = 0, kappahat = ones(25), randclone = randclone_1)
A = Population(lineages = [A1, A2], theta = 1/25)

V1 = Lineage(name = "virus", n = 1000, s = -1, kappa = ones(25), kappahat = ones(0))
V = Population(lineages = [V1], theta = 1/25)

evolve!(A, V, untilfixation = true )

# evolve two populations
A1 = Lineage(name = "antibody lin 1", n = 1000, s = 1, kappa = ones(25),
            shat = 1, kappahat = ones(25), randclone = randclone_1)
A = Population(lineages = [A1], theta = 1/25)

V1 = Lineage(name = "virus", n = 1000, s = -1, kappa = ones(25), kappahat = ones(0))
V = Population(lineages = [V1], theta = 1/25)

@time output = evolve!(A, V, tmax = 100)
