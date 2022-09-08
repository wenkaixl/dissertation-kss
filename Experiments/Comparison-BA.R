### PACKAGES
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
set.seed(19092022)

vertices = 20
edges_per_step = 1
n_sim_h0 = 500
n_sim_h1 = 500
B.list = c(200)
t_fun_selection = "comnb"

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

print("Barabasi")
print("Kernels:")
print(kernel.names)
print("B:")
print(B.list)
print("n:")
print(n_sim_h0)
print("t-function:")
print(t_fun_selection)

### EXPERIMENT
coef.h0 = 1
coef.h1.list = seq(0, 2, 0.25)

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

#t_fun based on common_neighbours lookup-table
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
  S=matrix(X,byrow=TRUE)[sample.index]
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

perform.test.sample = function(g.sim.data, generate.method, kernel, level = NULL,
                               coef.h1, n_sim_h1, B, sim.h0){
  
  N = length(sim.h0)
  pval.list = matrix(0,n_sim_h1)
  test.stats.list = matrix(0,n_sim_h1)
  d = dim(g.sim.data[[1]][,])[1]
  idx = sample.int(vertices^2, size = B, replace = TRUE)
  for (i in 1:n_sim_h1){
    X = g.sim.data[[i]][,]
    ##method
    test.out=GKSS.sampled(t_fun, X, sample.index = idx, kernel = kernel, level = level)
    test.stats = test.out[['stats.value']];
    pval = mean(rep(test.stats, times = N)<sim.h0)
    pval.list[i] = pval
    test.stats.list[i] = test.stats
  }
  list(pvalue=pval.list, stats=test.stats.list)
}

#simulate the null distribution
sim_gkss <- array(data = NA, dim = c(length(kernels), length(B.list), n_sim_h0))
g.sim <- vector("list", length = n_sim_h0)

for (i in 1:n_sim_h0) {
  G = sample_pa(n = vertices, m = edges_per_step, power = coef.h0, directed = FALSE)
  G = asNetwork(G)
  g.sim[[i]] <- G
}

### Calculate lookup-table from H0-simulations
if (t_fun_selection == "bideg") {
  t_fun = t_fun_bideg

  bideg.distribution.matrix = matrix(0, nrow = vertices, ncol = vertices)
  bideg.edge.distribution.matrix = matrix(0, nrow = vertices, ncol = vertices)

  for (G in g.sim) {
    G = asIgraph(G)
    x = G[,]

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

  for (G in g.sim) {
    G = asIgraph(G)
    x = G[,]

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
  for (G in g.sim) {
    G = asIgraph(G)

    density_sum = density_sum + edge_density(G, loops = FALSE)
  }
  density = density_sum / length(g.sim)
}

for (j in 1:length(B.list)){
  B = B.list[j]
  idx = sample.int(vertices^2, size = B, replace = TRUE) #sample indices
  for (i in 1:n_sim_h0) {
    G <- g.sim[[i]]
    X.h0 = G[,]
    for (k in 1:length(kernels)) {
      kernel = kernels[[k]]
      method.out=GKSS.sampled(t_fun, X.h0, sample.index = idx, kernel = kernel[[1]],
                              level = kernel[[2]])
      sim_gkss[k, j, i] = method.out[['stats.value']]
    }
    if(i%%10==0)print(paste("Iter",i, (sim_gkss[k,j,i])))
  }
}

l = length(B.list)
l.h1 = length(coef.h1.list)
test.stats.total = array(0, c(length(kernels), l ,l.h1, n_sim_h1))
power = array(0, c(length(kernels), l, l.h1))

powertest.parallel <- function(coef_idx, kernel_idx, B_idx) {
  coef = coef.h1.list[coef_idx]
  kernel = kernels[[kernel_idx]]
  B = B.list[B_idx]
  
  cat("Coefficient:", coef , "\n", sep = " ")
  
  coef.h1 = coef
  g.sim.h1 <- lapply(rep(vertices, n_sim_h1), sample_pa, m = edges_per_step, power = coef.h1, directed = FALSE)
  g.sim.h1 <- lapply(g.sim.h1, asNetwork)
  res = perform.test.sample(g.sim.h1, GKSS.sampled, kernel = kernel[[1]], level = kernel[[2]], coef.h1, n_sim_h1, B, sim_gkss[kernel_idx, B_idx,])
  pvalue = res[['pvalue']]
  test.stats.total[kernel_idx, B_idx, coef_idx, ] = res[['stats']]
  cat("Coefficient:", coef , "- finished\n",sep = " ")
  mean(pvalue < 0.025 | pvalue > 0.975)
}

for (kernel_idx in (1:length(kernels))) {
  for (B_idx in (1:length(B.list))) {
    B = B.list[B_idx]
    cat("--- Kernel:", kernel_idx, "B:", B, "\n", sep = " ")
    rejection_rates <- mclapply(as.list(1:l.h1), mc.cores = l.h1, FUN = powertest.parallel, kernel_idx=kernel_idx, B_idx=B_idx)
    power[kernel_idx, B_idx, ] <- unlist(rejection_rates)
    cat("--- Kernel:", kernel_idx, "B:", B, "- finished\n",sep = " ")
  }
}

dimnames(power) <- list(kernel.names, B.list, coef.h1.list)
df <- as.data.frame.table(power)
colnames(df) <- c("Kernel", "B", "Coefficient", "Power")
df[, 2] <- as.numeric(as.character(df[, 2]))
df[, 3] <- as.numeric(as.character(df[, 3]))
df[, 4] <- as.numeric(as.character(df[, 4]))
powerplot <- ggplot(data=df) + aes(x = Coefficient, y = Power, group = Kernel, colour = Kernel)+ geom_point() + geom_line() + ggtitle("Power of various kernels in BA graph (true coeff. = 1)")
print(df)

setwd("..")
setwd("Data")
#setwd("Barabasi-Albert Experiment")
ggsave("BAcomnb_all_1edges_v20_B200_n500.jpg", plot = powerplot)
write.csv(df,"BAcomnb_all_1edges_v20_B200_n500.csv", row.names = FALSE)