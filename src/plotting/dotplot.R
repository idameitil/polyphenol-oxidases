library(tidyverse)
library(ggdendro)
library(cowplot)
library(ggtree)
library(patchwork) 
library(ape)
library(stats)
library(grid)
library(gridExtra)

setwd("/Users/idamei/polyphenol-oxidases")

gene_cluster <- read_tsv('https://github.com/davemcg/davemcg.github.io/raw/master/content/post/scRNA_dotplot_data.tsv.gz')

markers <- gene_cluster$Gene %>% unique()

gene_cluster %>% filter(Gene %in% markers) %>% 
  mutate(`% Expressing` = (cell_exp_ct/cell_ct) * 100) %>% 
  ggplot(aes(x=cluster, y = Gene, color = count, size = `% Expressing`)) + 
  geom_point() 

# Mine 
df = read.table('data/mrbayes/all/clades/clades.csv', sep = ',', header=TRUE)

df = df %>%
  pivot_longer(!species, names_to = "group", values_to = "count")
df

species_unique <- df$species %>% unique()

gene_cluster %>% filter(Gene %in% markers) %>% 
  mutate(`% Expressing` = (cell_exp_ct/cell_ct) * 100) %>% 
  filter(count > 0, `% Expressing` > 1) %>% 
  ggplot(aes(x=cluster, y = Gene, color = count, size = `% Expressing`)) + 
  geom_point() + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank()) +
  scale_color_gradientn(colours = viridis::viridis(20), limits = c(0,4), oob = scales::squish, name = 'log2 (count + 1)')

tree <- read.tree('species.nwk')
#hc <- as.hclust.phylo(tree)
#dend <- as.dendrogram(hc)
dend <- chronos(tree)
ggtree_plot <- ggtree(dend) +
  geom_treescale() +
  geom_tiplab(align=TRUE, size=3) +
  xlim(0,1.1)
ggtree_plot <- ggtree(dend)

p2b = ggplot_build(ggtree_plot)

df = df %>% 
  mutate(species = factor(species, levels=p2b$data[[4]] %>% arrange(y) %>% pull(label)))

dotplot <- df %>% filter(species %in% species_unique) %>%
  filter(count > 0) %>% 
  ggplot(aes(x=group, y = species, color = group, size = count)) + 
  geom_point(show.legend = FALSE) +
  cowplot::theme_cowplot() + 
  geom_text(aes(label = count), color = 'white', size = 2) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.text.y = element_text(size = 6))

plot_grid(ggtree_plot, dotplot, nrow = 1, rel_widths = c(0.5,2), align = 'h')
grid.arrange(ggtree_plot, dotplot, ncol=2)
