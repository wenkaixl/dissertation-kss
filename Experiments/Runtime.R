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
vertices = 40
n_sim_h0 = 100
B.list = c(200)

#coef.h0 = c(-2,0) #sparse
coef.h0 = c(1,0) #dense

kernels = list(
  list(constant.kernel.sampled, NULL),
  list(wl.kernel.sampled, 1),
  list(wl.kernel.sampled, 3),
  list(wl.kernel.sampled, 5),
  list(sp.kernel.sampled, 1),
  list(grw.kernel.sampled, 1e-5),
  list(krw.kernel.sampled, 3),
  list(krw.kernel.sampled, 5),
  list(glet.kernel.sampled, 3),
  list(conglet.kernel.sampled, 3),
  list(conglet.kernel.sampled, 4),
  list(veh.kernel.sampled, 1)
)

kernel.names = c("Const.", "WL 1", "WL 3", "WL 5", "SP", "GRW", "KRW 3", "KRW 5", "GLET 3", "CONGLET 3", "CONGLET 4", "GVEH")

t_fun = function(x){
  Del = coef.h0[1] + coef.h0[2]*del_t2s(x)
  S = sigmoid(Del)
}

# construct the network  
un=network(vertices, directed = FALSE)
model0 <- un ~ edges + kstar(2)

# functions
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

library(microbenchmark)

#simulate the null distribution
g.sim <- simulate(model0, nsim=n_sim_h0,
                  coef=coef.h0,control=control.simulate(MCMC.burnin=1000+10*vertices, MCMC.interval=vertices))

time_overview <- array(data = NA, dim = c(length(B.list), length(kernels), n_sim_h0))
                       
for (i in (1:n_sim_h0)) {
 X.h0 = g.sim[[i]][,]
 for (j in (1:length(B.list))){
   B = B.list[j]
   idx = sample.int(vertices^2, size = B, replace = TRUE) #sample indices
   kernel_times = c()
   for (k in 1:length(kernels)) {
     kernel = kernels[[k]]
     time <- summary(microbenchmark(GKSS.sampled(t_fun, X.h0, sample.index = idx, kernel = kernel[[1]],
                             level = kernel[[2]]), times = 10, unit = "ms"))$median
     kernel_times[k] <- time
   }
   time_overview[j, , i] <- kernel_times
   if (i %% 1 == 0) print(paste("Iter", i))
 }
}

time_overview <- time_overview[1, , ]
time_overview
setwd("Evaluation time")
save(time_overview, file="time_sparse_n40_B200_sim100.rda")

dimnames(time_overview)[[1]] <- kernel.names

time_df <- as.data.frame.table(time_overview)
colnames(time_df) <- c("Kernel", "Iteration", "Time")
timeplot <- ggplot(data=time_df) + aes(x = Kernel, y = Time)+ geom_boxplot() + ggtitle("Time for calculation of sampled AgraSSt for B=200 in sparse E2S") + labs(y= "Time in ms", x = "Kernel")

setwd("..")
setwd("Data")
write.csv(time_df,"Const_time_dense_n40_B200_sim100_10_repeats.csv", row.names = FALSE)
ggsave("Const_Plot_time_dense_n40_B200_sim100_10_repeats.jpg", plot = timeplot)