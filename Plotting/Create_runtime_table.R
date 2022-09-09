library(xtable)

setwd("..")
setwd("Data")
df_s20 <- read.csv("./time_sparse_n20_B200_sim100.csv")
df_s40 <- read.csv("./time_sparse_n40_B200_sim100.csv")
df_d20 <- read.csv("./time_dense_n20_B200_sim100.csv")
df_d40 <- read.csv("./time_dense_n40_B200_sim100.csv")

colnames(df_complete) <- c("Kernel")

l <- list(df_s20, df_s40, df_d20, df_d40)
df_l <- list()

for (i in (1:4)) {
  df <- l[[i]]
  
  df_min <- aggregate(df[, c(1,3)], list(df[, 1]), FUN =min)
  df_min <- df_min[, c(1,3)]
  colnames(df_min) <- c("Kernel", "Min. runtime")
  df_mean <- aggregate(df[, c(1,3)], list(df[, 1]), FUN =mean)
  df_mean <- df_mean[, c(1,3)]
  colnames(df_mean) <- c("Kernel", "Avg. runtime")
  df_max <- aggregate(df[, c(1,3)], list(df[, 1]), FUN =max)
  df_max <- df_max[, c(1,3)]
  colnames(df_max) <- c("Kernel", "Max. runtime")
  df_merge = merge(df_min, df_mean, by = "Kernel")
  df_merge = merge(df_merge, df_max, by = "Kernel")
  
  df_l[[i]] <- df_merge
}

df_complete <- merge(df_l[[1]], df_l[[2]], by = "Kernel")
df_complete <- merge(df_complete, df_l[[3]], by = "Kernel")
df_complete <- merge(df_complete, df_l[[4]], by = "Kernel")
df_complete <- df_complete[order(df_complete$"Avg. runtime.y"),]

df_complete <- df_complete[,-(2:7)]

latex_table <- xtable(df_complete, digits = 2)
latex_table
