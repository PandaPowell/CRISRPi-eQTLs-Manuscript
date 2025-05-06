library(ggthemes)
library(ggplot2)

# Economist palette (light theme)
econ_colors <- ggthemes::ggthemes_data[["economist"]][["fg"]]
econ_dark = ggthemes::ggthemes_data[["economist"]][["bg"]]

# Convert to data frame for plotting
df <- bind_rows(econ_dark,econ_colors)

# Plot the colors
ggplot(df, aes(x = name, y = 1, fill = value)) +
  geom_tile() +
  geom_text(aes(label = name), angle = 90, vjust = -0.5, hjust = 0) +
  scale_fill_identity() +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  ggtitle("Economist Color Palette (Light Theme)")
