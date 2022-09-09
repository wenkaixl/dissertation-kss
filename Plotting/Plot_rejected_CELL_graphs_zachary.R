library(igraph)

vertices = 34

get_graph <- function(vec) {
  m <- matrix(vec, nrow = vertices)
  G <- graph_from_adjacency_matrix(m, mode = "undirected")
  G
}

df_zach <- read.csv(file.choose(), header = FALSE)
df_sim <- read.csv(file.choose())
df_rej_idx <- read.csv(file.choose())

rej_idx <- as.logical(df_rej_idx[df_rej_idx[,1] == "GVEH 0.1", 2:101]) & as.logical(df_rej_idx[df_rej_idx[,1] == "GVEH 1", 2:101])
idx <- which(rej_idx)
nrej_idx <- as.logical((1-df_rej_idx[1, 2:101]) & (1-df_rej_idx[2, 2:101]) & (1-df_rej_idx[3, 2:101]) & (1-df_rej_idx[4, 2:101]) & (1-df_rej_idx[5, 2:101]) & (1-df_rej_idx[6, 2:101]) & (1-df_rej_idx[7, 2:101]) & (1-df_rej_idx[8, 2:101]) & (1-df_rej_idx[9, 2:101]))
non_rej_idx <- which(nrej_idx)

G_train <- get_graph(df_zach[1,])
zach_layout <- layout_(G_train, nicely())
zach_deg <- degree(G_train, mode="all")

cfg <- cluster_fast_greedy(G_train)
cfg <- cutat(cfg, 2)
vertex_attr(G_train, "Faction") <- cfg


par(mfrow = c(3,3))
plot(G_train, vertex.color = c("red", "pink")[cfg], vertex.label = NA, vertex.size = 15, edge.width = 3, edge.color = "black",
     layout = zach_layout)

for (i in (59*100+1):(59*100+8)) {
  G_sim <- get_graph(df_sim[i,])
  sim_deg <- degree(G_sim, mode = "all")
  sim_cfg <- cluster_fast_greedy(G_sim)
  sim_cfg <- cutat(sim_cfg, 2)
  if (sim_cfg[33] == 1) sim_cfg <- 3 - sim_cfg
  plot(G_sim, vertex.color = c("blue", "lightblue")[sim_cfg], vertex.label = NA, vertex.size = 15, edge.width = 3, edge.color = "black",
       layout = zach_layout)
}
