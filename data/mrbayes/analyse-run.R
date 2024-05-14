library(ggplot2)

setwd("/Users/idamei/polyphenol-oxidases/data/mrbayes/all")
df = read.table(file = 'out.nex.mcmc', sep = '\t', header = TRUE, skip=5)
df['GenMio'] = df['Gen'] / 1000000
ggplot(df, aes(x=GenMio, y=AvgStdDev.s.)) +
  geom_point()

setwd("/Users/idamei/hpc/fungal-10-mixed")
df = read.table(file = 'test.nex.mcmc', sep = '\t', header = TRUE, skip=5)
df['GenMio'] = df['Gen'] / 1000000
ggplot(df, aes(x=GenMio, y=AvgStdDev.s.)) +
  geom_point() +
  ylim(0,0.2)

setwd("/Users/idamei/fungal-10-mixed")
df = read.table(file = 'test.nex.mcmc', sep = '\t', header = TRUE, skip=5)
df['GenMio'] = df['Gen'] / 1000000
ggplot(df, aes(x=GenMio, y=AvgStdDev.s.)) +
  geom_point() +
  ylim(0,0.2)
