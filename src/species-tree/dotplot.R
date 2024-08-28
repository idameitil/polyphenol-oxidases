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
df = read.table('data/epa-ng/filtered-out/clades.csv', sep = ',', header=TRUE)
df2 = read.table('data/epa-ng/fungi-with-lignin-degraders/clades.csv', sep = ',', header=TRUE)

names(df)[names(df) == 's'] = 'u'
names(df2)[names(df2) == 's'] = 'u'

species_unique <- df$species %>% unique()
species_unique2 <- df2$species %>% unique()

df = df %>%
  pivot_longer(!species, names_to = "group", values_to = "count")
df2 = df2 %>%
  pivot_longer(!species, names_to = "group", values_to = "count")

# Read metadata
metadata <- read_excel('data/proteome-tree/class-representatives.xlsx')
metadata$species <- gsub(' ', '_', metadata$species)

# Read metadata2
metadata2 <- read_excel('data/proteome-tree/fungal-order-representatives-andlignindegraders.xlsx')
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

#'#C87EE9', '#A2F8F9'
colors_groups_all = c('a'='#6985b5', 'b'='#8e1730', 'c'='#52099b', 'd'='#ecd75a', 'e'='#F056EA', 'f'='#fb8b34', 'g'='#00057C', 'h'='#119d58', 'i'='#10aefd', 'j'='#db9758', 'k'='#116F6D', 'l'='#f71252', 'u'='#000000')

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
  geom_point(show.legend = FALSE, alpha=2) +
  geom_text(aes(label = count), color = 'white', size = 1.5) +
  theme_minimal() +
  theme(
        axis.text.x = element_text(colour = colors_groups_all, size = 14),
        axis.text.y = element_text(size = 7, hjust=0, face = 'italic'),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        plot.margin = unit(c(0.5,0.5,0.5,3), 'cm')
        ) +
  scale_color_manual(values = colors_groups_all) +
  scale_y_discrete(labels = change_species_name) +
  scale_size_continuous(range= c(1,6)) +
  coord_cartesian(clip = "off")  # Prevent clipping
plot_grid(ggtree_plot, dotplot, nrow = 1, rel_widths = c(2.7,3), align = 'h')

colors_groups = colors_groups = c('d'='#ecd75a','e'= '#F056EA',  'f'='#fb8b34', 'k'='#116F6D', 'l'='#f71252', 'u'='#000000')
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
  scale_y_discrete(labels = change_species_name, limits=rev) +
  coord_cartesian(clip = "off")  # Prevent clipping
dotplot2

# Save pdf
pdf(file="manuscript/figures/dotplot.pdf",
    width = 8,
    height = 11)
plot_grid(ggtree_plot, dotplot, nrow = 1, rel_widths = c(3,4.5), align = 'h')
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
