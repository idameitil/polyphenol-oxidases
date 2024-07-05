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
#df = read.table('data/species-tree/clades2.csv', sep = ',', header=TRUE)
df = read.table('data/mrbayes/all-seeds-0619/clades/clades.csv', sep = ',', header=TRUE)

names(df)[names(df) == 'singletons'] = 's'

df = df %>%
  pivot_longer(!species, names_to = "group", values_to = "count")

species_unique <- df$species %>% unique()

# Read metadata
metadata <- read_excel('data/proteome-tree/class-representatives.xlsx')
metadata$species <- gsub(' ', '_', metadata$species)

# Merge
merged_table <- left_join(df, metadata, by = "species")

# Read tree
tree <- read.tree('data/species-tree/species.nwk')

# Use tree order to order dotplot
tree_plot_labelled <- ggtree(tree) +
  geom_tiplab()
p2b = ggplot_build(tree_plot_labelled)
species_order = p2b$data[[3]] %>% arrange(y) %>% pull(label)
merged_table = merged_table %>% 
  mutate(species = factor(species, levels=species_order))

# Plot tree
ggtree_plot <- ggtree(tree) +
  theme(plot.margin = unit(c(0.5,0,0.5,0), 'cm'))

#'#C87EE9', '#A2F8F9', '#00057C', '#116F6D', '#ecd75a'
colors_groups = c('#10aefd', '#db9758', '#f71252', '#119d58', '#6985b5', '#8e1730', '#52099b', '#fb8b34', '#F056EA', '#000000')

change_species_name <- function(x) gsub("(^\\w)\\w*_(\\w+)", "\\1. \\2", x)
# Make dotplot
dotplot <- merged_table %>% filter(species %in% species_unique) %>%
  filter(count > 0) %>% 
  ggplot(mapping = aes(
    x=group, 
    y = species, 
    color = group, 
    size = count)) + 
  geom_point(show.legend = FALSE, alpha=1) +
  geom_text(aes(label = count), color = 'white', size = 1.5) +
  theme_minimal() +
  theme(
        axis.text.x = element_text(colour = colors_groups, size = 14),
        axis.text.y = element_text(size = 7, hjust=0, face = 'italic'),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        plot.margin = unit(c(0.5,0.5,0.5,2), 'cm')
        ) +
  scale_color_manual(values = colors_groups) +
  scale_y_discrete(labels = change_species_name)
plot_grid(ggtree_plot, dotplot, nrow = 1, rel_widths = c(4,3), align = 'h')

# Save pdf
pdf(file="manuscript/figures/dotplot.pdf",
    width = 9,
    height = 8)
plot_grid(ggtree_plot, dotplot, nrow = 1, rel_widths = c(4,3), align = 'h')
dev.off()

# Write phylum list
phylum_ordered <- species_order %>%
  as.data.frame() %>%
  rename(species = ".") %>%
  left_join(metadata, by = "species") %>%
  pull(phylum)

file_name <- file("data/species-tree/phylum.txt")
writeLines(phylum_ordered, file_name)
close(file_name)

