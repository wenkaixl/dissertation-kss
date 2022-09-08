library(ergm)
library(network)
library(igraph)
library(graphkernels)

transposed.index = function(w, vertices){
  col_idx <- ( (w-1) %% vertices) + 1
  row_idx <- ((w-1) %/% vertices) + 1
  transposed_index <- (col_idx - 1) * vertices + row_idx
  transposed_index
}

# creates a list of graphs where exactly one edge indicator is flipped
compute.transition.list = function(X){
  P=list()
  for (w in 1:length(X)){
    x = X
    x[w] = abs(1-X[w])
    w_t = transposed.index(w, nrow(X))
    if (w_t != w) x[w_t] = abs(1-X[w])
    G = graph_from_adjacency_matrix(x, mode = "undirected")
    P[[w]] = G
  }
  P[[length(X)+1]] = graph_from_adjacency_matrix(X, mode = "undirected")
  P
}

compute.sampled.list = function(X, sample.index){
  P=list()
  l = length(sample.index)
  for (i in 1:l){
    x = X
    w = sample.index[i]
    x[w] = abs(1 - X[w])
    w_t = transposed.index(w, nrow(X))
    if (w_t != w) x[w_t] = abs(1-X[w])
    G = graph_from_adjacency_matrix(x, mode = "undirected")
    P[[w]] = G
  }
  P[[l+1]] = graph_from_adjacency_matrix(X, mode = "undirected")
  P
}

compute.normalized = function(K){
  V = diag(K)
  D.mat = diag(1.0/sqrt(V))
  D.mat%*%K%*%D.mat
}

## Weisfeiler_Lehman Graph Kernel
wl.kernel.sampled=function(X, sample.index, level=3, diag=1, normalize=TRUE){
  n = length(sample.index)
  P = compute.sampled.list(X, sample.index)
  kernel.matrix = CalculateWLKernel(P, level)
  rm(P) #removes P (takes a lot of memory)
  K = kernel.matrix[1:n,1:n] + kernel.matrix[n+1,n+1] #k(x-i flipped, x-j flipped) + k(x,x)
  K.vec = kernel.matrix[1:n,n+1] #kernel between flipped and original graph
  K = K + outer(K.vec, K.vec, function(x,y)x+y) #outer= k(x-i flipped, x) + k(x, x-j flipped)
  #K(i,j) = k(x,x)+k(x-i flipped, x)+k(x,x-j flipped)+k(x-i flipped, x-j flipped)
  if(normalize)K = compute.normalized(K)
  if(diag==0)diag(K)<-0
  K
}

## k-Step Random Walk Kernel
krw.kernel.sampled=function(X, sample.index, level = 3, diag=1, normalize=TRUE){
  n = length(sample.index)
  P = compute.sampled.list(X, sample.index)
  kernel.matrix = CalculateKStepRandomWalkKernel(P, level)
  rm(P) 
  K = kernel.matrix[1:n,1:n] + kernel.matrix[n+1,n+1]
  K.vec = kernel.matrix[1:n,n+1]
  K = K + outer(K.vec, K.vec, function(x,y)x+y)
  if(normalize)K = compute.normalized(K)
  if(diag==0)diag(K)<-0
  K
}

## Geometric Random Walk Kernel
grw.kernel.sampled=function(X, sample.index, level=3, diag=1, normalize=TRUE){
  n = length(sample.index)
  P = compute.sampled.list(X, sample.index)
  kernel.matrix = CalculateGeometricRandomWalkKernel(P, level)
  rm(P) 
  K = kernel.matrix[1:n,1:n] + kernel.matrix[n+1,n+1]
  K.vec = kernel.matrix[1:n,n+1]
  K = K + outer(K.vec, K.vec, function(x,y)x+y)
  if(normalize)K = compute.normalized(K)
  if(diag==0)diag(K)<-0
  K
}

## Shortest Path Kernel
sp.kernel.sampled=function(X, sample.index, level=1, diag=1, normalize=TRUE){
  n = length(sample.index)
  P = compute.sampled.list(X, sample.index)
  kernel.matrix = CalculateShortestPathKernel(P)
  rm(P)
  K = kernel.matrix[1:n,1:n] + kernel.matrix[n+1,n+1]
  K.vec = kernel.matrix[1:n,n+1]
  K = K + outer(K.vec, K.vec, function(x,y)x+y)
  if(normalize)K = compute.normalized(K)
  if(diag==0)diag(K)<-0
  K
}

## Generalised Shortest Path Kernel - utility functions
## Generalised Shortest Path Kernel
compute.gsp.matrix=function(g) {
  vertices = length(V(g))
  
  sp.length <- shortest.paths(g)
  sp.count <- matrix(NA, nrow = vertices, ncol = vertices)
  
  for (i in (1:vertices)) {
    sp.info <- all_shortest_paths(g, i)
    sp.count[i,] <- sp.info$nrgeo
  }
  list(sp.length, sp.count)
}

CalculateGeneralisedShortestPathKernel=function(P) {
  n = length(P)
  K <- matrix(NA, nrow = n, ncol = n)
  for (i in (1:n)) {
    m1 <- compute.gsp.matrix(P[[i]])
    length1 <- m1[[1]]
    count1 <- m1[[2]]
    for (j in (i:n)) {
      m2 <- compute.gsp.matrix(P[[j]])
      length2 <- m2[[1]]
      count2 <- m2[[2]]
      
      comb <- length1 * length2 * (count1 == count2)
      K[i,j] <- sum(comb)
      K[j,i] <- K[i,j]
    }
  }
  K
}

## Generalised Shortest Path Kernel
gsp.kernel.sampled=function(X, sample.index, level=1, diag=1, normalize=TRUE){
  n = length(sample.index)
  P = compute.sampled.list(X, sample.index)
  kernel.matrix = CalculateGeneralisedShortestPathKernel(P)
  rm(P)
  K = kernel.matrix[1:n,1:n] + kernel.matrix[n+1,n+1]
  K.vec = kernel.matrix[1:n,n+1]
  K = K + outer(K.vec, K.vec, function(x,y)x+y)
  if(normalize)K = compute.normalized(K)
  if(diag==0)diag(K)<-0
  K
}

## Graphlet Kernel
glet.kernel.sampled=function(X, sample.index, level=3, diag=1, normalize=TRUE){
  n = length(sample.index)
  P = compute.sampled.list(X, sample.index)
  kernel.matrix = CalculateGraphletKernel(P, level)
  rm(P) 
  K = kernel.matrix[1:n,1:n] + kernel.matrix[n+1,n+1]
  K.vec = kernel.matrix[1:n,n+1]
  K = K + outer(K.vec, K.vec, function(x,y)x+y)
  if(normalize)K = compute.normalized(K)
  if(diag==0)diag(K)<-0
  K
}

## Connected Graphlet Kernel
conglet.kernel.sampled=function(X, sample.index, level=3, diag=1, normalize=TRUE){
  n = length(sample.index)
  P = compute.sampled.list(X, sample.index)
  kernel.matrix = CalculateConnectedGraphletKernel(P, level)
  rm(P) 
  K = kernel.matrix[1:n,1:n] + kernel.matrix[n+1,n+1]
  K.vec = kernel.matrix[1:n,n+1]
  K = K + outer(K.vec, K.vec, function(x,y)x+y)
  if(normalize)K = compute.normalized(K)
  if(diag==0)diag(K)<-0
  K
}

## Vertex-Edge Histogram Kernel
veh.kernel.sampled=function(X, sample.index, level=0.1, diag=1, normalize=TRUE){
  n = length(sample.index)
  P = compute.sampled.list(X, sample.index)

  for (i in (1:length(P))) {
    g = P[[i]]
    P[[i]] = set_vertex_attr(g, "degree", index = V(g), value=degree(g, v=V(g), loops = FALSE))
  }

  kernel.matrix = CalculateVertexEdgeHistGaussKernel(P,level)
  rm(P) 
  K = kernel.matrix[1:n,1:n] + kernel.matrix[n+1,n+1]
  K.vec = kernel.matrix[1:n,n+1]
  K = K + outer(K.vec, K.vec, function(x,y)x+y)
  if(normalize)K = compute.normalized(K)
  if(diag==0)diag(K)<-0
  K
}

## Constant Kernel
constant.kernel.sampled=function(X, sample.index, diag=1, normalize = TRUE){
  n = length(sample.index)
  K = matrix(0,ncol=n,nrow = n) + 1
  if(diag==0)diag(K)<-0
  K
}

## Constant Kernel x1000 (only for sanity check)
constant.kernel.sampled1000=function(X, sample.index, diag=1, normalize = TRUE){
  n = length(sample.index)
  K = matrix(0,ncol=n,nrow = n) + 1000
  if(diag==0)diag(K)<-0
  K
}

## Gaussian Twostars Kernel
# This is a clunkly implementation as there are faster ways to count two-stars if only one edge flips
# Kernel not used in dissertation

CalculateTwoStarsKernel=function(P, lambda = 1) {
  n = length(P)
  twostars <- rep(NA, n)
  for (i in 1:n) {
    G <- as_adjacency_matrix(P[[i]], sparse = FALSE)
    twostar_matrix <- G %*% t(G)
    diag(twostar_matrix) <- 0
    vertices = dim(twostar_matrix)[1]
    W = matrix(1, nrow = vertices, ncol = 1)
    twostars[i] <- t(W)%*%twostar_matrix%*%W / 2
  }
  twostars_diff_sq <- outer(twostars, twostars, function(x,y)(x-y)^2)
  kernel.matrix <- exp(-lambda*(1/2)*twostars_diff_sq)
  kernel.matrix
}

gaussian.twostars.kernel.sampled=function(X, sample.index, diag=1, normalize=TRUE){
  n = length(sample.index)
  P = compute.sampled.list(X, sample.index)
  kernel.matrix = CalculateTwoStarsKernel(P)
  rm(P) 
  K = kernel.matrix[1:n,1:n] + kernel.matrix[n+1,n+1]
  K.vec = kernel.matrix[1:n,n+1]
  K = K + outer(K.vec, K.vec, function(x,y)x+y)
  if(normalize)K = compute.normalized(K)
  if(diag==0)diag(K)<-0
  K
}