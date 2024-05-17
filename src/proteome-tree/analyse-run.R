library(ggplot2)
library(ggpubr)

par(mfrow=c(2,1))

setwd("/Users/idamei/polyphenol-oxidases/data/mrbayes/all/run")
df = read.table(file = 'out.nex.mcmc', sep = '\t', header = TRUE, skip=5)
df['GenMio'] = df['Gen'] / 1000000
all <- ggplot(df, aes(x=GenMio, y=AvgStdDev.s.)) +
  geom_point() +
  ylim(0,0.2)

setwd("/Users/idamei/hpc/fungal-10-mixed")
df = read.table(file = 'test.nex.mcmc', sep = '\t', header = TRUE, skip=5)
df['GenMio'] = df['Gen'] / 1000000
fungal <- ggplot(df, aes(x=GenMio, y=AvgStdDev.s.)) +
  geom_point() +
  ylim(0,0.2)

setwd("/Users/idamei/hpc/fungal-new")
df = read.table(file = 'test.nex.mcmc', sep = '\t', header = TRUE, skip=5)
df['GenMio'] = df['Gen'] / 1000000
fungal_new <- ggplot(df, aes(x=GenMio, y=AvgStdDev.s.)) +
  geom_point() +
  ylim(0,0.2)

ggarrange(all, fungal, fungal_new,
          labels = c("A", "B", "C"),
          ncol = 3, nrow = 1)
