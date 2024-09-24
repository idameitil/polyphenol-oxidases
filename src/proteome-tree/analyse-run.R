library(ggplot2)
library(ggpubr)

setwd("/Users/idamei/polyphenol-oxidases/data/mrbayes/all/run")
df = read.table(file = 'out.nex.mcmc', sep = '\t', header = TRUE, skip=5)
df['GenMio'] = df['Gen'] / 1000000
all <- ggplot(df, aes(x=GenMio, y=AvgStdDev.s.)) +
  geom_point() +
  ylim(0,0.2) +
  xlim(0, 60)


setwd("/Users/idamei/hpc/all-seeds-0619")
df = read.table(file = 'out.nex.mcmc', sep = '\t', header = TRUE, skip=5)
df['GenMio'] = df['Gen'] / 1000000
all_seeds <- ggplot(df, aes(x=GenMio, y=AvgStdDev.s.)) +
  geom_point() +
  ylim(0,0.2) +
  xlim(0, 90)

setwd("/Users/idamei/hpc/short-fungal-0507")
df = read.table(file = 'out.nex.mcmc', sep = '\t', header = TRUE, skip=5)
df['GenMio'] = df['Gen'] / 1000000
short_fungal <- ggplot(df, aes(x=GenMio, y=AvgStdDev.s.)) +
  geom_point() +
  ylim(0,0.2) +
  xlim(0, 70)

setwd("/Users/idamei/0816")
df = read.table(file = 'out.nex.mcmc', sep = '\t', header = TRUE, skip=5)
df['GenMio'] = df['Gen'] / 1000000
all <- ggplot(df, aes(x=GenMio, y=AvgStdDev.s.)) +
  geom_point() +
  ylim(0,0.2) +
  xlim(0, 210)

setwd("/Users/idamei/polyphenol-oxidases/")
df = read.table(file = 'out.nex.mcmc', sep = '\t', header = TRUE, skip=5)
df['GenMio'] = df['Gen'] / 1000000
hemo <- ggplot(df, aes(x=GenMio, y=AvgStdDev.s.)) +
  geom_point() +
  ylim(0,0.2) +
  xlim(0, 330)

setwd("/Users/idamei/polyphenol-oxidases/strucaligned")
df = read.table(file = 'out.nex.mcmc', sep = '\t', header = TRUE, skip=5)
df['GenMio'] = df['Gen'] / 1000000
struc <- ggplot(df, aes(x=GenMio, y=AvgStdDev.s.)) +
  geom_point() +
  ylim(0,0.3) +
  xlim(0, 20)

ggarrange(hemo, struc,
          labels = c("A", "B"),
          ncol = 2, nrow = 1)

