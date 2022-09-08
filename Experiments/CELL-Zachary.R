## PACKAGES
library(ergm)
library(network)
library(igraph)
library(graphkernels)
library(ggplot2)
library(parallel)
library(intergraph)

setwd("..")
setwd("Utilities")
source("kernels_thesis.R")
source("Base_fun.R")

### PARAMETERS FOR EXPERIMENT
vertices = 34
n_gen = 100
n_rep = 100
B.list = c(200)
t_fun_selection = "density"

set.seed(19082022)

kernels = list(
  list(constant.kernel.sampled, NULL),
  list(wl.kernel.sampled, 1),
  list(wl.kernel.sampled, 3),
  list(wl.kernel.sampled, 5),
  list(sp.kernel.sampled, 1),
  list(grw.kernel.sampled, 1e-5),
  list(grw.kernel.sampled, 1e-4),
  list(grw.kernel.sampled, 1e-3),
  list(grw.kernel.sampled, 1e-2),
  list(grw.kernel.sampled, 5e-2),
  list(krw.kernel.sampled, 2),
  list(krw.kernel.sampled, 3),
  list(krw.kernel.sampled, 4),
  list(krw.kernel.sampled, 5),
  list(glet.kernel.sampled, 3),
  list(conglet.kernel.sampled, 3),
  list(conglet.kernel.sampled, 4),
  list(veh.kernel.sampled, 0.1),
  list(veh.kernel.sampled, 1),
  list(veh.kernel.sampled, 10),
  list(veh.kernel.sampled, 100)
)

kernel.names = c("Const.", "WL 1", "WL 3", "WL 5", "SP", "GRW 1e-5", "GRW 1e-4", "GRW 1e-3", "GRW 1e-2", "GRW 5e-2",
 "KRW 2", "KRW 3", "KRW 4", "KRW 5", "GLET 3", "CONGLET 3", "CONGLET 4",
 "GVEH 0.1", "GVEH 1", "GVEH 10", "GVEH 100")

print("Zachary CELL")
print("Kernels:")
print(kernel.names)
print("B:")
print(B.list)
print("n:")
print(n_gen)
print("t-function:")
print(t_fun_selection)

#load data
setwd("..")
setwd("Data")
zach_adj_vec = read.csv("Zachary_adjacency.csv", header = FALSE)
zach_adj = matrix(zach_adj_vec[1,], nrow = vertices)
zach = graph_from_adjacency_matrix(zach_adj, mode = "undirected")

zach.sim = read.csv("CELL_graphs_from_zachary.csv", header = FALSE)

### EXPERIMENT

#t_fun based on (deg(i), deg(j))-lookup table
t_fun_bideg = function(x){
  G = graph_from_adjacency_matrix(x, mode = "undirected")
  deg = degree(G, loops = FALSE)
  density = edge_density(G, loops = FALSE)

  q.matrix = matrix(NA, nrow = vertices, ncol = vertices)
  for (i in (1:vertices)){
    for (j in (1:vertices)){
      deg_i = deg[i]
      deg_j = deg[j]

      #if there exists edge between i and j, then their degree is one lower in the graph, where
      #edge indicator (i,j) is missing
      if (x[i,j] == 1){
        deg_i = deg_i - 1
        deg_j = deg_j - 1
      }

      if (bideg.distribution.matrix[(deg_i + 1), (deg_j + 1)] == 0) {
        q.matrix[i,j] = density
      } else {
        q.matrix[i,j] = (bideg.edge.distribution.matrix[(deg_i + 1), (deg_j + 1)]) / (bideg.distribution.matrix[(deg_i + 1), (deg_j + 1)])
      }
    }
  }
  q.matrix
}

t_fun_density = function(x) {
  q.matrix = matrix(density, nrow = vertices, ncol = vertices)
  q.matrix
}

#t_fun based on #common_neighbours lookup-table
t_fun_comnb = function(x){
  G = graph_from_adjacency_matrix(x, mode = "undirected")
  density = edge_density(G, loops = FALSE)
  q.matrix = matrix(NA, nrow = vertices, ncol = vertices)
  for (i in (1:vertices)){
    for (j in (1:vertices)){
      comnb = length(setdiff(intersect(nb_list[[i]], nb_list[[j]]), c(i,j)))

      if (comnb.distribution.matrix[comnb + 1] == 0) {
        q.matrix[i,j] = density #no other pair of vertices with the same number of common neighbours is observed
      } else {
        q.matrix[i,j] = (comnb.edge.distribution.matrix[comnb + 1]) / (comnb.distribution.matrix[comnb + 1])
      }
    }
  }
  q.matrix
}

# functions
GKSS.sampled=function(t_fun, X, sample.index, kernel, level = NULL,
                      normalize=TRUE, v.scale=TRUE, diag=1)
{
  S = matrix(X,byrow=TRUE)[sample.index]
  n.sim=length(S)
  S.t =  t_fun(X)[sample.index]
  S.t.vec=abs(S - S.t)
  S.mat = S.t.vec %*% t(S.t.vec)
  Ky.mat = (S*2-1)%*%t(S*2-1)
  
  n = length(sample.index)
  
  if (is.null(level)) {
    K = kernel(X, sample.index, diag = diag, normalize = normalize)
  } else {
    K = kernel(X, sample.index, level = level, diag = diag, normalize = normalize)
  }
  
  J.kernel = S.mat * Ky.mat * K
  
  W=rep(1/n,n)
  J.kernel.out=J.kernel
  if(diag==0)diag(J.kernel)=rep(0,n.sim)
  v.kernel = var(S)
  stats.value=t(W)%*%J.kernel%*%W 
  if(v.scale)nKSD = stats.value* sqrt(v.kernel)
  #Return:
  #stats.value: n times KSD
  #J.kernel: J kernel matrix for wild bootstrapt 
  list(stats.value=stats.value,J.kernel=J.kernel.out, K=K, S=S.mat)
}

perform.test.cell.sample = function(g.sim, generate.method, kernel, level = NULL,
                               model.h1, coef.h1, n_sim_h1, B, sim.h0){
  
  gen.gkss.list = matrix(0,n_gen)
  for (i in (1:n_gen)){
    idx = sample.int(vertices^2, size = B, replace = TRUE)
    X = g.sim[i,]
    X = apply(X, 2, as.numeric)
    X = matrix(X, byrow = TRUE, nrow = vertices)
    X = graph_from_adjacency_matrix(X, mode = "undirected")
    X = asNetwork(X)
    X = X[,]
    ##method
    test.out=GKSS.sampled(t_fun, X, sample.index = idx, kernel = kernel, level = level)
    test.stats = test.out[['stats.value']]
    gen.gkss.list[i] = test.stats
  }
  pval = mean(gen.gkss.list < sim.h0)
  reject = (pval < 0.025 | pval > 0.975)
  list(pval,reject, gen.gkss.list)
}

### Calculate lookup-table from H0-simulations
if (t_fun_selection == "bideg") {
  t_fun = t_fun_bideg

  bideg.distribution.matrix = matrix(0, nrow = vertices, ncol = vertices)
  bideg.edge.distribution.matrix = matrix(0, nrow = vertices, ncol = vertices)

  for (i in (1:(n_rep*n_gen))) {
    x = zach.sim[i, ]
    x = matrix(x, nrow = vertices, ncol = vertices, byrow = TRUE)
    G = graph_from_adjacency_matrix(x, mode = "undirected")

    deg = degree(G, loops = FALSE)
    for (i in (1:(vertices-1))){
      for (j in ((i+1):vertices)){
        deg_i = deg[i]
        deg_j = deg[j]

        #if there exists edge between i and j, then their degree is one lower in the graph, where
        #edge indicator (i,j) is missing
        if(x[i,j] == 1) {
          deg_i = deg_i - 1
          deg_j = deg_j - 1
        }

        bideg.distribution.matrix[(deg_i+1), (deg_j+1)] = bideg.distribution.matrix[(deg_i+1), (deg_j+1)] + 1
        if (deg_i != deg_j) bideg.distribution.matrix[(deg_j+1), (deg_i+1)] = bideg.distribution.matrix[(deg_j+1), (deg_i+1)] + 1
        if (x[i,j] == 1) {
          bideg.edge.distribution.matrix[(deg_i+1), (deg_j+1)] = bideg.edge.distribution.matrix[(deg_i+1), (deg_j+1)] + 1
          if (deg_i != deg_j) bideg.edge.distribution.matrix[(deg_j+1), (deg_i+1)] = bideg.edge.distribution.matrix[(deg_j+1), (deg_i+1)] + 1
        }
      }
    }
  }
}

if (t_fun_selection == "comnb") {
  t_fun = t_fun_comnb

  comnb.distribution.matrix = matrix(0, nrow = 1, ncol = vertices)
  comnb.edge.distribution.matrix = matrix(0, nrow = 1, ncol = vertices)

  for (i in (1:(n_rep*n_gen))) {
    x = zach.sim[i, ]
    x = matrix(x, nrow = vertices, ncol = vertices, byrow = TRUE)
    G = graph_from_adjacency_matrix(x, mode = "undirected")

    nb_list = neighborhood(G, order = 1, nodes = V(G))
    for (i in (1:(vertices-1))){
      for (j in ((i+1):vertices)){
        comnb = length(setdiff(intersect(nb_list[[i]], nb_list[[j]]), c(i,j)))

        comnb.distribution.matrix[comnb + 1] = comnb.distribution.matrix[comnb + 1] + 1 
        if (x[i,j] == 1) {
          comnb.edge.distribution.matrix[comnb + 1] = comnb.edge.distribution.matrix[comnb + 1] + 1
        }
      }
    }
  }
}

if (t_fun_selection == "density") {
  t_fun = t_fun_density

  density_sum = 0
  for (i in (1:(n_rep*n_gen))) {
    x = zach.sim[i, ]
    x = matrix(x, nrow = vertices, ncol = vertices, byrow = TRUE)
    G = graph_from_adjacency_matrix(x, mode = "undirected")

    density_sum = density_sum + edge_density(G, loops = FALSE)
  }
  density = density_sum / (n_rep*n_gen)
}

#simulate the null distribution
sim_gkss <- array(data = NA, dim = c(length(kernels), length(B.list), n_rep))
g.sim <- zach

for (j in 1:length(B.list)){
  B = B.list[j]
  g.sim = asNetwork(g.sim)
  X = g.sim[,]
  X = apply(X, 2, as.numeric)
  X = matrix(X, byrow = TRUE, nrow = vertices)
  X = graph_from_adjacency_matrix(X, mode = "undirected")
  X = asNetwork(X)
  X = X[,]
  X.h0 = X
  for (k in 1:length(kernels)) {
    kernel = kernels[[k]]
    for (i in (1:n_rep)) {
      idx = sample.int(vertices^2, size = B, replace = TRUE) #sample indices
      method.out=GKSS.sampled(t_fun, X.h0, sample.index = idx, kernel = kernel[[1]],
                              level = kernel[[2]])
      sim_gkss[k, j, i] = method.out[['stats.value']]
    }
  }
}

l = length(B.list)
power = array(0, c(length(kernels), l))
rejected_indices_df = data.frame(matrix(NA, nrow = length(kernels), ncol = n_rep)) #if B.list only contains one value

powertest.parallel.cell <- function(rep_nr, kernel_idx, B_idx) {
  kernel = kernels[[kernel_idx]]
  B = B.list[B_idx]
  
  if (rep_nr %% 10 == 0) cat("Repetition:", rep_nr , "\n", sep = " ")
  
  g.sim <- zach.sim[(((rep_nr-1)*n_gen + 1):(rep_nr*n_gen)), ]
  res = perform.test.cell.sample(g.sim, GKSS.sampled, kernel = kernel[[1]], level = kernel[[2]], model0, coef.h1, n_sim_h1, B, sim_gkss[kernel_idx, B_idx, rep_nr])
  if (rep_nr %% 10 == 0) cat("Repetition:", rep_nr , "- finished\n",sep = " ")
  #if (rep_nr==1) print(res[[3]])
  #if (rep_nr==1) print(sim_gkss)
  res[[2]]
}

for (kernel_idx in (1:length(kernels))) {
  for (B_idx in (1:length(B.list))) {
    B = B.list[B_idx]
    cat("--- Kernel:", kernel_idx, "B:", B, "\n", sep = " ")
    rejections <- mclapply(as.list(1:n_rep), mc.cores = 16, powertest.parallel.cell, kernel_idx=kernel_idx, B_idx=B_idx)
    rejection_rate <- mean(unlist(rejections))
    #print(rejection_rates)
    power[kernel_idx, B_idx] <- rejection_rate
    rejected_indices_df[kernel_idx, ] <- unlist(rejections)
    #print(power)
    #print(class(power))
    cat("--- Kernel:", kernel_idx, "B:", B, "- finished\n",sep = " ")
  }
}

dimnames(power) <- list(kernel.names, B.list)
df <- as.data.frame.table(power)
colnames(df) <- c("Kernel", "B", "Power")
powerplot <- ggplot(data=df) + aes(x = Kernel, y = Power)+ geom_col() + ggtitle("Rejection probability of kernels on CELL-generated graphs for B = 200 for Karate club network")
powerplot <- powerplot + ylab("Rejection probability for two-sided test at 5%-level")

print(df)

rownames(rejected_indices_df) <- kernel.names
colnames(rejected_indices_df) <- (1:n_gen)

#setwd("CELL Zachary experiment")
ggsave("n100_density_all_CELL_zach_power-all.jpg", plot = powerplot)
write.csv(df,"n100_density_all_CELL_zach_power-all.csv", row.names = FALSE)
write.csv(rejected_indices_df,"n100_density_all_CELL_zach_rejected_batches-all.csv")