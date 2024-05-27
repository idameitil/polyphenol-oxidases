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
library("readxl")
library(dplyr)
library(ggtext)

setwd("/Users/idamei/polyphenol-oxidases")

# Read clades table
df = read.table('data/species-tree/clades2.csv', sep = ',', header=TRUE)

names(df)[names(df) == 'singletons'] = 's'

df = df %>%
  pivot_longer(!species, names_to = "group", values_to = "count")
df

species_unique <- df$species %>% unique()

# Read metadata
metadata <- read_excel('data/proteome-tree/class-representatives.xlsx')
metadata$species <- gsub(' ', '_', metadata$species)

# Merge
merged_table <- left_join(df, metadata, by = "species")

# Tree
## Method 1
#x <- read.dendrogram('data/species-tree/species.nwk')
#p1 <- plot(x, horiz = TRUE)

## Method 2
tree <- read.tree('data/species-tree/species.nwk')
dend <- chronos(tree)
ggtree_plot <- ggtree(dend) +
  geom_treescale() +
  geom_tiplab(align=TRUE, size=3) +
  xlim(0,1.1)

## Use tree order to order dotplot
p2b = ggplot_build(ggtree_plot)
species_order = p2b$data[[6]] %>% arrange(y) %>% pull(label)
merged_table = merged_table %>% 
  mutate(species = factor(species, levels=species_order))

## Method 3
ggtree_plot <- ggtree(tree) +
  theme(plot.margin = unit(c(0.5,0,0.5,0.5), 'cm')) 

#'#C87EE9'
#'#A2F8F9'
#'#00057C'
#'#116F6D'
colors_groups = c('#10aefd', '#db9758', '#f71252', '#119d58', '#6985b5', '#ecd75a', '#8e1730', '#52099b', '#fb8b34', '#F056EA', '#000000')

myfun <- function(x) gsub("(^\\w)\\w*_(\\w+)", "\\1. \\2", x)
myfun('Homo_sapiens')
# Make dotplot
dotplot <- merged_table %>% filter(species %in% species_unique) %>%
  filter(count > 0) %>% 
  ggplot(mapping = aes(
    x=group, 
    y = species, 
    color = group, 
    size = count)) + 
  geom_point(show.legend = FALSE, alpha=1) +
#  cowplot::theme_cowplot() + 
  geom_text(aes(label = count), color = 'white', size = 1.5) +
  theme_minimal() +
  theme(
        #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.x = element_text(colour = colors_groups, size = 14),
        axis.text.y = element_text(size = 9, hjust=0, face = 'italic'),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        plot.margin = unit(c(0.5,0.5,0.5,3.3), 'cm')
        ) +
  scale_color_manual(values = colors_groups) +
  scale_y_discrete(labels = myfun)
plot_grid(ggtree_plot, dotplot, nrow = 1, rel_widths = c(1,3), align = 'h')

# Arange
## Method 1
pdf(file="manuscript/figures/dotplot.pdf",
    width = 6,
    height = 8)
plot_grid(ggtree_plot, dotplot, nrow = 1, rel_widths = c(1,3), align = 'h')
dev.off()
## Method 2
#grid.arrange(ggtree_plot, dotplot, ncol=2)




# Phylum stuff
df_filtered <- merged_table %>% filter(species %in% species_unique) %>%
  filter(count > 0)

phylum_ordered <- species_order %>%
  as.data.frame() %>%
  rename(species = ".") %>%
  left_join(metadata, by = "species") %>%
  pull(phylum)
#output for python
#cat(phylum_ordered, sep="', '")

result_table <- data.frame(
  phylum = phylum_ordered,
  x = seq(length(phylum_ordered)),
  y = factor(seq(length(phylum_ordered)))
)

result_table %>%
  ggplot(aes(x, y, color=phylum)) +
  geom_blank() +
  theme(axis.text.x = element_text(size=6)) +
  scale_y_discrete(labels = phylum_ordered) +
  theme(plot.margin = margin(0, 0, 0, 0))

+
  theme(axis.text.y = element_text(size=6),
        axis.line.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank())

