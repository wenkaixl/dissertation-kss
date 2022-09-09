library(ggplot2)

kernel.names <- c("Const.", "GRW 1e-5", "GRW 1e-4", "GRW 1e-3", "GRW 1e-2", "GRW 5e-2")
two_rows <- TRUE
glue <- FALSE
glue.kernel.names <- c("Const.")

setwd("..")
setwd("Data")
df <- read.csv(file.choose())
df <- subset(df, Kernel %in% kernel.names)
df$Kernel <-factor(df$Kernel, levels=kernel.names)
if (glue == TRUE) {
  df2 <- read.csv(file.choose())
  df2 <- subset(df2, Kernel %in% glue.kernel.names)
  df <- rbind(df2, df)
  df$Kernel <- factor(df$Kernel, levels=c(glue.kernel.names, kernel.names))
}
powerplot <- ggplot(data=df) + aes(x = Coefficient, y = Power, group = Kernel, colour = Kernel)
powerplot <- powerplot + ylim(0,1)
powerplot <- powerplot + scale_x_continuous(limits=c(0.1, 0.5), breaks = c(0.1, 0.3, 0.5))
powerplot <- powerplot + geom_point(size = 3)
powerplot <- powerplot + geom_line(linetype = "dashed", size = 1.3)
powerplot <- powerplot + ylab("Rejection rate")
powerplot <- powerplot + xlab("Parameter of alternative model")
powerplot <- powerplot + theme(legend.position="bottom")
powerplot <- powerplot + theme(text = element_text(size=25))
if (two_rows == TRUE)powerplot <- powerplot + guides(colour=guide_legend(nrow=2,byrow=FALSE))

powerplot
#ggsave("Dense_GRW_n20_B200.jpg", plot = powerplot)