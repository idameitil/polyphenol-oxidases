library(ggplot2)
library(ggpubr)

setwd("/Users/idamei/polyphenol-oxidases/data/mrbayes/all/run")
df = read.table(file = 'out.nex.mcmc', sep = '\t', header = TRUE, skip=5)
df['GenMio'] = df['Gen'] / 1000000
all <- ggplot(df, aes(x=GenMio, y=AvgStdDev.s.)) +
  geom_point() +
  ylim(0,0.2) +
  xlim(0, 60)

setwd("/Users/idamei/hpc/fungal-10-mixed")
df = read.table(file = 'test.nex.mcmc', sep = '\t', header = TRUE, skip=5)
df['GenMio'] = df['Gen'] / 1000000
fungal <- ggplot(df, aes(x=GenMio, y=AvgStdDev.s.)) +
  geom_point() +
  ylim(0,0.2) +
  xlim(0, 60)

setwd("/Users/idamei/hpc/fungal-newnew")
df = read.table(file = 'test.nex.mcmc', sep = '\t', header = TRUE, skip=5)
df['GenMio'] = df['Gen'] / 1000000
fungal_new <- ggplot(df, aes(x=GenMio, y=AvgStdDev.s.)) +
  geom_point() +
  ylim(0,0.2) +
  xlim(0, 75)

setwd("/Users/idamei/hpc/all-new")
df = read.table(file = 'out.nex.mcmc', sep = '\t', header = TRUE, skip=5)
df['GenMio'] = df['Gen'] / 1000000
all_new <- ggplot(df, aes(x=GenMio, y=AvgStdDev.s.)) +
  geom_point() +
  ylim(0,0.2) +
  xlim(0, 70)

setwd("/Users/idamei/hpc/all-seeds-0619")
df = read.table(file = 'out.nex.mcmc', sep = '\t', header = TRUE, skip=5)
df['GenMio'] = df['Gen'] / 1000000
all_seeds <- ggplot(df, aes(x=GenMio, y=AvgStdDev.s.)) +
  geom_point() +
  ylim(0,0.2) +
  xlim(0, 70)

ggarrange(all, all_new, all_seeds
          labels = c("A", "B", "C"),
          ncol = 3, nrow = 1)
