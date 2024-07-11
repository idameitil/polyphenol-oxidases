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
df2 = read.table('data/fungi-genome-figure/clades.csv', sep = ',', header=TRUE)

names(df)[names(df) == 's'] = 'u'
names(df2)[names(df2) == 's'] = 'u'

df = df %>%
  pivot_longer(!species, names_to = "group", values_to = "count")
df2 = df2 %>%
  pivot_longer(!species, names_to = "group", values_to = "count")

species_unique <- df$species %>% unique()
species_unique2 <- df2$species %>% unique()

# Read metadata
metadata <- read_excel('data/proteome-tree/class-representatives.xlsx')
metadata$species <- gsub(' ', '_', metadata$species)

# Read metadata2
metadata2 <- read_excel('data/proteome-tree/fungal-order-representatives.xlsx')
metadata2$species <- gsub(' ', '_', metadata2$species)

# Merge
merged_table <- left_join(df, metadata, by = "species")
merged_table2 <- left_join(df2, metadata2, by = "species")

# Read tree
tree <- read.tree('data/species-tree/species.nwk')

# Use tree order to order dotplot
tree_plot_labelled <- ggtree(tree) +
  geom_tiplab()
p2b = ggplot_build(tree_plot_labelled)
species_order = p2b$data[[3]] %>% arrange(y) %>% pull(label)
merged_table = merged_table %>% 
  mutate(species = factor(species, levels=species_order))
merged_table2 = merged_table2 %>% 
  mutate(species = factor(species, levels=metadata2$species))

# Plot tree
ggtree_plot <- ggtree(tree) +
  theme(plot.margin = unit(c(0.5,0,0.5,0), 'cm'))

#'#C87EE9', '#A2F8F9', '#00057C', '#116F6D', '#ecd75a'
colors_groups = c('#10aefd', '#db9758', '#f71252', '#119d58', '#6985b5', '#8e1730', '#52099b', '#fb8b34', '#F056EA', '#000000')

#change_species_name <- function(x) gsub("(^\\w)\\w*_(\\w+)", "\\1. \\2", x)
change_species_name <- function(x) gsub("(^\\w+)\\w*_(\\w+)", "\\1 \\2", x)
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
plot_grid(ggtree_plot, dotplot, nrow = 1, rel_widths = c(3,3), align = 'h')

colors_groups = c('#f71252', '#fb8b34', '#F056EA', '#000000')
dotplot2 <- merged_table2 %>% filter(species %in% species_unique2) %>%
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
plot_grid(ggtree_plot, dotplot2, nrow = 1, rel_widths = c(3,3), align = 'h')

# Save pdf
pdf(file="manuscript/figures/dotplot.pdf",
    width = 8,
    height = 8)
plot_grid(ggtree_plot, dotplot, nrow = 1, rel_widths = c(3,3.5), align = 'h')
dev.off()

# Save pdf 2
pdf(file="manuscript/figures/dotplot2.pdf",
    width = 4,
    height = 8)
plot_grid(dotplot2, nrow = 1, rel_widths = c(3,3), align = 'h')
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

# Write phylum list 2
wide_table <- merged_table2 %>% 
  pivot_wider(names_from = group, values_from  = count) %>% 
  arrange(factor(species, levels = metadata2$species))
wide_table$phylum

file_name <- file("data/fungi-genome-figure/phylum.txt")
writeLines(wide_table$phylum, file_name)
close(file_name)

file_name <- file("data/fungi-genome-figure/class.txt")
writeLines(wide_table$class, file_name)
close(file_name)

file_name <- file("data/fungi-genome-figure/order.txt")
writeLines(wide_table$order, file_name)
close(file_name)

file_name <- file("data/fungi-genome-figure/family.txt")
writeLines(wide_table$family, file_name)
close(file_name)

file_name <- file("data/fungi-genome-figure/genus.txt")
writeLines(wide_table$genus, file_name)
close(file_name)
