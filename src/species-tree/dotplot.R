library(tidyverse)
library(ggdendro)
library(cowplot)
library(ggtree)
library(patchwork) 
library(ape)
library(stats)
library(grid)
library(gridExtra)
library("phylogram")

setwd("/Users/idamei/polyphenol-oxidases")

df = read.table('clades2.csv', sep = ',', header=TRUE)

df = df %>%
  pivot_longer(!species, names_to = "group", values_to = "count")
df

species_unique <- df$species %>% unique()

# Tree
## Method 1
x <- read.dendrogram('data/species-tree/species.nwk')
p1 <- plot(x, horiz = TRUE)

## Method 2
tree <- read.tree('species.nwk')
dend <- chronos(tree)
ggtree_plot <- ggtree(dend) +
  geom_treescale() +
  geom_tiplab(align=TRUE, size=3) +
  xlim(0,1.1)

## Use tree order to order dotplot
p2b = ggplot_build(ggtree_plot)
df = df %>% 
  mutate(species = factor(species, levels=p2b$data[[6]] %>% arrange(y) %>% pull(label)))

## Method 3
ggtree_plot <- ggtree(dend)

#C87EE9
#A2F8F9
#00057C
colors_groups = c('#10aefd', '#db9758', '#f71252', '#119d58', '#6985b5', '#ecd75a', '#fb8b34', '#8e1730', '#52099b', '#F056EA', '#116F6D')

# Make dotplot
dotplot <- df %>% filter(species %in% species_unique) %>%
  filter(count > 0) %>% 
  ggplot(aes(x=group, y = species, color = group, size = count)) + 
  geom_point(show.legend = FALSE) +
#  cowplot::theme_cowplot() + 
  geom_text(aes(label = count), color = 'white', size = 1.5) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size = 6),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank()) +
  scale_color_manual(values = colors_groups)

# Arange
## Method 1
plot_grid(ggtree_plot, dotplot, nrow = 1, rel_widths = c(2,2), align = 'h')
## Method 2
grid.arrange(ggtree_plot, dotplot, ncol=2)
