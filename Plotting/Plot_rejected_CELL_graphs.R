library(igraph)

vertices = 20

get_graph <- function(vec) {
  m <- matrix(vec, nrow = vertices)
  G <- graph_from_adjacency_matrix(m, mode = "undirected")
  G
}

df_train <- read.csv(file.choose())
df_sim <- read.csv(file.choose())
df_rej_idx <- read.csv(file.choose())

rej_idx <- as.logical(df_rej_idx[df_rej_idx[,1] == "GVEH 0.1", 2:101])
idx <- which(rej_idx)

for (i in idx) {
  G_train <- get_graph(df_train[i,])
  G_sim <- get_graph(df_sim[i,])
  par(mfrow = c(1,2))
  layout_train = layout_(G_sim, nicely())
  #layout_train = layout_(G_train, nicely())
  plot(G_train, vertex.color = "red", vertex.label = NA, vertex.size = 12, edge.width = 3, edge.color = "black",
       layout = layout_train)
  plot(G_sim, vertex.color = "blue", vertex.label = NA, vertex.size = 12, edge.width = 3, edge.color = "black",
       layout = layout_train)
  readline(prompt="Press [enter] to continue")
}