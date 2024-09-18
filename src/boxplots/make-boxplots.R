library(dplyr)
library(ggplot2)
library(readr)

setwd("/Users/idamei/polyphenol-oxidases")

data <- read_delim("data/boxplots/boxplot-data.tsv", na = c("", "NA"), delim='\t')

na_df = data[is.na(data$class),]

data$class[data$order == 'Mycoplasmoidales'] <- 'Mycoplasmoidales'
data$class[data$order == 'Polyangiales'] <- 'Polyangiales'
data$class[data$order == 'Peronosporales'] <- 'Peronosporales'
data$class[data$order == 'Saprolegniales'] <- 'Saprolegniales'
data$class[data$order == 'Polyangiales'] <- 'Polyangiales'
data$class[data$order == 'Testudines'] <- 'Testudines'

unique_values_in_order <- unique(data$class[match(data$class, data$class)])
data$class <- factor(data$class, levels = unique_values_in_order)

data$count_tyrosinases[is.na(data$count_tyrosinases)] <- 0

# Count the number of members in each class
class_counts <- data %>%
  group_by(class) %>%
  summarise(count = n())

# Filter out classes with fewer than 5 members
filtered_classes <- class_counts %>%
  filter(count >= 5) %>%
  select('class')

filtered_data <- data %>%
  filter(class %in% filtered_classes$class & !is.na(class) & class != '')

# Step 1: Count the number of members per class and include whether the 'kingdom' column has NA
# class_counts <- data %>%
#   group_by(class) %>%
#   summarise(
#     count = n(),
#     kingdom_na_count = sum(is.na(kingdom))
#   )
# 
# # Step 2: Filter out classes with fewer than 5 members only if 'kingdom' is NA
# filtered_classes <- class_counts %>%
#   filter(count >= 10 | kingdom_na_count == 0) %>%
#   select('class')
# 
# # Step 3: Filter the original data to keep only the classes that met the criteria
# filtered_data <- data %>%
#   filter(class %in% filtered_classes$class)

y_breaks <- seq(0, max(filtered_data$count_tyrosinases, na.rm = TRUE), by = 5)

# Create a boxplot with the filtered data
p1 <- ggplot(filtered_data, aes(x = class, y = count_tyrosinases)) +
  geom_boxplot() +
  geom_jitter(size = 0.1, width = 0.3, height = 0, alpha = 0.7, color = "blue") +
  xlab('Taxonomic class') +
  ylab('Number of PPO genes') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_y_continuous(breaks = y_breaks)
p1

pdf('manuscript/figures/boxplot.pdf', width = 20, height = 11)
p1
dev.off()


# flipped
p1 <- ggplot(filtered_data, aes(x = class, y = count_tyrosinases)) +
  geom_boxplot(color = 'grey50') +
  geom_jitter(size = 0.1, width = 0.3, height = 0.1, alpha = 0.7, color = "blue") +
  xlab('Taxonomic class') +
  ylab('Number of PPO genes') +
  theme(axis.text.x = element_text()) +
  scale_y_continuous(breaks = y_breaks) +
  coord_flip()
p1

filtered_data <- filtered_data %>%
  mutate(Taxonomy = case_when(
    kingdom == 'Fungi' ~ 'Fungi',
    kingdom == 'Oomycota' ~ 'Oomycota',
    phylum == 'Chordata' ~ 'Chordata',
    phylum == 'Cnidaria' ~ 'Cnidaria',
    phylum %in% c('Bacteroidota', 'Planctomycetota', 'Pseudomonadota', 'Acidobacteriota',
                  'Chloroflexota', 'Cyanobacteriota', 'Bacillota', 'Actinomycetota', 
                    'Deinococcota', 'Gemmatimonadota', 'Armatinomadetes', 
                  'Verrucomicrobiota', 'Nitrospirota', 'Aquificota', 'Balneolota',
                  'Bdellovibrionota', 'Campylobacterota', 'Chlamydiota', 'Chlorobiota',
                  'Deferribacterota', 'Fusobacteriota', 'Ignavibacteriota', 'Mycoplasmatota',
                  'Myxococcota', 'Rhodothermota', 'Spirochaetota', 'Synergistota', 'Thermodesulfobacteriota',
                  'Thermotogota') ~ 'Bacteria',
    class %in% c('Deltaproteobacteria') ~ 'Bacteria',
    phylum == 'Mollusca' ~ 'Mollusc',
    phylum == 'Oomycota' ~ 'Oomycota',
    kingdom == 'Viridiplantae' ~ 'Plant',
    phylum == 'Brachiopoda' ~ 'Brachiopod',
    phylum == 'Nematoda' ~ 'Nematode',
    phylum == 'Arthropoda' ~ 'Arthropod',
    phylum %in% c('Thermoproteota', 'Euryarchaeota', 'Korarchaeota', 'Nanoarchaeota', 
                  'Thaumarchaeota', 'Candidatus Thermoplasmatota', 'Nitrososphaerota') ~ 'Archaea',
    TRUE ~ 'Other'  # Default case for any other values
  ))

# Load required libraries
library(ggplot2)
library(dplyr)

filtered_data$Taxonomy <- factor(filtered_data$Taxonomy, levels = c('Chordata', 'Mollusc', 'Arthropod', 
                                                                    'Nematode', 'Cnidaria', 'Fungi',
                                                                    'Oomycota', 'Plant', 'Bacteria',
                                                                    'Archaea', 'Other'))

sorted_data <- filtered_data %>%
  arrange(Taxonomy, phylum)

# Step 2: Update the factor levels based on the sorted order
sorted_data <- sorted_data %>%
  mutate(class = factor(class, levels = unique(class)))

levels(sorted_data$class)

# Define the custom colors as a named vector (R equivalent of a dictionary)
value2color <- c(
  'Plant' = '#00FF00',
  'Other' = '#FFFFFF',
  'Fungi' = '#964B00',
  'Chordata' = '#ffff00',
  'Mollusc' = '#FF00FF',
  'Cnidaria' = '#ff0000',
  'Oomycota' = '#937584',
  'Bacteria' = '#00ffff',
  'Nematode' = '#008080',
  'Arthropod' = '#FF7F50',
  'Archaea' = '#0000FF'
)

# Create the plot with custom colors
p1 <- ggplot(sorted_data, aes(x = class, y = count_tyrosinases)) +
  geom_boxplot(color = 'grey50') +
  geom_jitter(size = 0.1, width = 0.3, height = 0.1, alpha = 0.7, color = "blue") +
  xlab('Taxonomic class') +
  ylab('Number of PPO genes') +
  theme(axis.text.x = element_text()) +
  scale_y_continuous(limits = c(-1.5, 85), breaks = y_breaks, expand = c(0, 0)) +  # Remove extra padding
  # Add color strip based on the 'adapted' column and use custom colors
  geom_tile(aes(x = class, y = -1, fill = Taxonomy), width = 1, height = 1) +
  scale_fill_manual(values = value2color) +  # Apply custom fill colors for the color strip
  scale_color_manual(values = value2color) +  # Apply custom colors for jitter points
  coord_flip(clip = "off") +  # Keep flipped and allow tiles outside plot
  theme(
    plot.margin = unit(c(1, 1, 2, 1), "lines")  # Adjust margins if needed
  ) +
  xlab('Taxonomic class') +
  ylab('Number of PPO genes')
p1
other <- filtered_data[filtered_data['Taxonomy'] == 'Other',]

dev.off()
pdf('manuscript/figures/boxplot_flipped.pdf', width = 11, height = 15)
p1
dev.off()
