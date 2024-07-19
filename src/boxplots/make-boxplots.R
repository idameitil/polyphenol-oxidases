library(dplyr)
library(ggplot2)

setwd("/Users/idamei/polyphenol-oxidases")

data = read.table('data/boxplots/boxplot-data.tsv', sep = '\t', header=TRUE)

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

# Filter the original data to keep only classes with 5 or more members
filtered_data <- data %>%
  filter(class %in% filtered_classes$class & !is.na(class) & class != '')

y_breaks <- seq(0, max(filtered_data$count_tyrosinases, na.rm = TRUE), by = 5)

# Create a boxplot with the filtered data
p1 <- ggplot(filtered_data, aes(x = class, y = count_tyrosinases)) +
  geom_boxplot() +
  geom_jitter(size = 0.1, width = 0.3, height = 0, alpha = 0.7, color = "blue") +
  xlab('Taxonomic class') +
  ylab('Number of PPO genes') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_continuous(breaks = y_breaks)

pdf('boxplot.pdf', width = 20, height = 11)
p1
dev.off()
