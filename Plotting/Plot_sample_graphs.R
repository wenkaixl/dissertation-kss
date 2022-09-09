library(igraph)
library(ergm)
library(network)
library(intergraph)

vertices = 20

# E2S
un=network(vertices, directed = FALSE)
model0 <- un ~ edges + kstar(2)

g1 <- simulate(model0, nsim=1, coef=c(-0,-0.5),control=control.simulate(MCMC.burnin=1000+10*vertices, MCMC.interval=vertices))
g2 <- simulate(model0, nsim=1, coef=c(-0,-0.2),control=control.simulate(MCMC.burnin=1000+10*vertices, MCMC.interval=vertices))
g3 <- simulate(model0, nsim=1, coef=c(-0,0),control=control.simulate(MCMC.burnin=1000+10*vertices, MCMC.interval=vertices))
g4 <- simulate(model0, nsim=1, coef=c(-0,0.5),control=control.simulate(MCMC.burnin=1000+10*vertices, MCMC.interval=vertices))

g1 <- asIgraph(g1)
g2 <- asIgraph(g2)
g3 <- asIgraph(g3)
g4 <- asIgraph(g4)

plot(g1, layout=layout_nicely, vertex.label = NA, vertex.color = "royalblue3", vertex.size = 7, edge.width = 1.3)
plot(g2, layout=layout_nicely, vertex.label = NA, vertex.color = "royalblue3", vertex.size = 7, edge.width = 1.3)
plot(g3, layout=layout_nicely, vertex.label = NA, vertex.color = "royalblue3", vertex.size = 7, edge.width = 1.3)
plot(g4, layout=layout_nicely, vertex.label = NA, vertex.color = "royalblue3", vertex.size = 7, edge.width = 1.3)

# GRG
grg1 <- sample_grg(nodes = vertices, radius = 0.1, torus = TRUE)
grg2 <- sample_grg(nodes = vertices, radius = 0.2, torus = TRUE)
grg3 <- sample_grg(nodes = vertices, radius = 0.3, torus = TRUE)
grg4 <- sample_grg(nodes = vertices, radius = 0.4, torus = TRUE)
grg5 <- sample_grg(nodes = vertices, radius = 0.5, torus = TRUE)

plot(grg1, layout=layout_nicely, vertex.label = NA, vertex.color = "royalblue3", vertex.size = 7, edge.width = 1.3)
plot(grg2, layout=layout_nicely, vertex.label = NA, vertex.color = "royalblue3", vertex.size = 7, edge.width = 1.3)
plot(grg3, layout=layout_nicely, vertex.label = NA, vertex.color = "royalblue3", vertex.size = 7, edge.width = 1.3)
plot(grg4, layout=layout_nicely, vertex.label = NA, vertex.color = "royalblue3", vertex.size = 7, edge.width = 1.3)
plot(grg5, layout=layout_nicely, vertex.label = NA, vertex.color = "royalblue3", vertex.size = 7, edge.width = 1.3)

# BA-1
ba1 <- sample_pa(n = vertices, power = 0, m = 1, directed = FALSE)
ba2 <- sample_pa(n = vertices, power = 0.5, m = 1, directed = FALSE)
ba3 <- sample_pa(n = vertices, power = 1, m = 1, directed = FALSE)
ba4 <- sample_pa(n = vertices, power = 1.5, m = 1, directed = FALSE)
ba5 <- sample_pa(n = vertices, power = 2, m = 1, directed = FALSE)

plot(ba1, layout=layout_nicely, vertex.label = NA, vertex.color = "royalblue3", vertex.size = 7, edge.width = 1.3)
plot(ba2, layout=layout_nicely, vertex.label = NA, vertex.color = "royalblue3", vertex.size = 7, edge.width = 1.3)
plot(ba3, layout=layout_nicely, vertex.label = NA, vertex.color = "royalblue3", vertex.size = 7, edge.width = 1.3)
plot(ba4, layout=layout_nicely, vertex.label = NA, vertex.color = "royalblue3", vertex.size = 7, edge.width = 1.3)
plot(ba5, layout=layout_nicely, vertex.label = NA, vertex.color = "royalblue3", vertex.size = 7, edge.width = 1.3)

# BA-2
ba2_1 <- sample_pa(n = vertices, power = 0, m = 2, directed = FALSE)
ba2_2 <- sample_pa(n = vertices, power = 0.5, m = 2, directed = FALSE)
ba2_3 <- sample_pa(n = vertices, power = 1, m = 2, directed = FALSE)
ba2_4 <- sample_pa(n = vertices, power = 1.5, m = 2, directed = FALSE)
ba2_5 <- sample_pa(n = vertices, power = 2, m = 2, directed = FALSE)

plot(ba2_1, layout=layout_nicely, vertex.label = NA, vertex.color = "royalblue3", vertex.size = 7, edge.width = 1.3)
plot(ba2_2, layout=layout_nicely, vertex.label = NA, vertex.color = "royalblue3", vertex.size = 7, edge.width = 1.3)
plot(ba2_3, layout=layout_nicely, vertex.label = NA, vertex.color = "royalblue3", vertex.size = 7, edge.width = 1.3)
plot(ba2_4, layout=layout_nicely, vertex.label = NA, vertex.color = "royalblue3", vertex.size = 7, edge.width = 1.3)
plot(ba2_5, layout=layout_nicely, vertex.label = NA, vertex.color = "royalblue3", vertex.size = 7, edge.width = 1.3)

